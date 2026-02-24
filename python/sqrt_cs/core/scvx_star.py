"""
SCvx* — Successive Convex Approximation solver.

Implements Algorithm 1 from the paper, with Augmented Lagrangian penalty
and physics-informed trust regions.

Reference:
  K. Oguri, "Successive Convex Approximation with Feasibility Guarantee
  via Augmented Lagrangian for Non-Convex Optimal Control Problems", 2023.
"""

from __future__ import annotations

import time
from typing import Any

import cvxpy as cp
import numpy as np

from sqrt_cs.core.scp_params import SCPParams
from sqrt_cs.core.scp_iter import SCPIter
from sqrt_cs.core.scp_problem import SCPProblem
from sqrt_cs.problems.sqrt_qr_cs import SqrtQRCovarianceSteering


class SCvxStar:
    """
    SCP solver for covariance steering problems.

    Usage
    -----
    solver = SCvxStar(problem, params)
    result = solver.solve()
    """

    def __init__(self, problem: SCPProblem, params: SCPParams | None = None):
        self.prob   = problem
        self.params = params or SCPParams()
        self.iter_state: SCPIter | None = None
        self.solve_time: float = 0.0

    # ================================================================== #
    # Main solve loop
    # ================================================================== #

    def solve(self) -> dict[str, np.ndarray]:
        """
        Run the SCP loop until convergence or k_max iterations.

        Returns the solution dict (same format as SqrtQRCovarianceSteering.extract_vals).
        """
        p  = self.params
        pr = self.prob

        # ---- Initial reference ----
        ref = pr.initial_ref()

        # ---- Evaluate initial residuals to size the slack vectors ----
        g_ref = pr.noncvx_eq(ref)
        h_ref = pr.noncvx_ineq(ref)
        J0_ref = self._eval_objective_at_ref(ref)

        n_eq   = len(g_ref)
        n_ineq = len(h_ref)

        # obstacle / inexact constraints sized dynamically each iter
        self.iter_state = SCPIter(
            n_eq=n_eq,
            n_ineq=n_ineq,
            n_inexact_ineq=0,
            r_init=p.r_init,
            w_init=p.w_init,
        )
        it = self.iter_state

        # ---- Define CVXPY variables once ----
        cvx_vars = pr.define_vars()

        if p.verbose:
            self._print_header()

        t_start = time.perf_counter()

        for k in range(1, p.k_max + 1):
            pr.pre_iteration(k, ref)

            # ---- Build and solve convex subproblem ----
            t_sub = time.perf_counter()
            status, vals, J0_k, slack_eq, slack_ineq, slack_inexact = \
                self._solve_subproblem(cvx_vars, ref, it)
            t_sub = time.perf_counter() - t_sub

            if vals is None:
                if p.verbose:
                    print(f"  iter {k:4d} | solver status: {status} — skipping step")
                # shrink trust region and try again
                it.r = max(it.r / p.alpha1, p.r_min)
                continue

            pr.post_iteration(k, vals, ref)

            # ---- Compute actual constraint violations at new point ----
            g_k = pr.noncvx_eq(vals)
            h_k = pr.noncvx_ineq(vals)
            chi_k = np.linalg.norm(np.concatenate([g_k, np.maximum(h_k, 0)]))

            # ---- Penalties ----
            P_ref = self._penalty(g_ref, h_ref, slack_eq * 0, slack_ineq * 0, it)
            P_k   = self._penalty(g_k,   h_k,   slack_eq,     slack_ineq,     it)

            delta_J = (J0_ref + P_ref) - (J0_k + P_k)    # actual improvement
            L_k = J0_k + float(np.sum(it.lam * slack_eq)) + \
                         float(np.sum(it.mu_al * np.maximum(slack_ineq, 0))) + \
                         (it.w / 2) * (float(np.sum(slack_eq**2)) +
                                       float(np.sum(np.maximum(slack_ineq, 0)**2)))
            delta_L = (J0_ref + P_ref) - L_k              # predicted improvement
            rho = delta_J / (delta_L + 1e-30)

            # ---- Accept / reject step ----
            ref_updated = False
            mul_updated = False

            if rho >= p.rho0:
                ref         = vals
                g_ref       = g_k
                h_ref       = h_k
                J0_ref      = J0_k
                ref_updated = True

                # Multiplier update criterion (Eq. 15 in SCvx*)
                if abs(delta_J) < min(it.delta, chi_k):
                    it.lam   = it.lam   + it.w * g_k
                    it.mu_al = np.maximum(it.mu_al + it.w * h_k, 0.0)
                    it.w     = min(p.beta * it.w, p.w_max)
                    it.delta = max(p.gamma_al * it.delta, p.tol_opt)
                    mul_updated = True

            # ---- Trust region update ----
            if rho < p.rho1:
                it.r = max(it.r / p.alpha1, p.r_min)
            elif rho >= p.rho2:
                it.r = min(it.r * p.alpha2, p.r_max)

            it.record(J0_k, chi_k, delta_J, delta_L, rho, ref_updated, mul_updated)

            if p.verbose:
                self._print_iter(k, J0_k, chi_k, delta_J, rho, it.r, it.w, ref_updated, mul_updated)

            # ---- Convergence check ----
            dz = self._variable_change_norm(vals, ref if ref_updated else ref)
            if (abs(delta_J) <= p.tol_opt
                    and chi_k <= p.tol_feas
                    and dz     <= p.tol_change):
                if p.verbose:
                    print(f"\n  Converged at iteration {k}.")
                break
        else:
            if p.verbose:
                print(f"\n  Reached k_max = {p.k_max}.")

        self.solve_time = time.perf_counter() - t_start

        # ---- Post-process solution ----
        pr.postprocess(ref, {})
        return ref

    # ================================================================== #
    # Convex subproblem
    # ================================================================== #

    def _solve_subproblem(
        self,
        cvx_vars: dict,
        ref: dict[str, np.ndarray],
        it: SCPIter,
    ) -> tuple[str, dict | None, float, np.ndarray, np.ndarray, np.ndarray]:
        """
        Build and solve one convex subproblem.

        Returns
        -------
        status      : solver status string
        vals        : extracted variable values  (None if infeasible)
        J0          : objective value
        slack_eq    : equality slack values
        slack_ineq  : inequality slack values
        slack_inexact: inexact ineq slack values
        """
        p  = self.params
        pr = self.prob

        # ---- Warm-start variables from reference ----
        self._warm_start(cvx_vars, ref)

        # ---- Slack variables for nonconvex constraints ----
        n_eq   = it.n_eq
        n_ineq = it.n_ineq

        sl_eq     = cp.Variable(n_eq,   name="slack_eq")
        sl_ineq   = cp.Variable(n_ineq, name="slack_ineq", nonneg=True) if n_ineq > 0 else None
        sl_inexact = None   # sized dynamically from convexified_inexact_ineq

        # ---- Constraint lists ----
        constraints: list[cp.Constraint] = []

        # Always-convex equality constraints
        constraints += pr.convex_eq(cvx_vars)

        # Always-convex inequality constraints
        constraints += pr.convex_ineq(cvx_vars)

        # Linearized nonconvex equalities  g(z) ≈ 0  →  g_lin(z) = slack_eq
        noncvx_eq_lin = pr.noncvx_eq_relaxed(cvx_vars, ref)
        if noncvx_eq_lin and n_eq > 0:
            # The relaxed constraints should equal the slack
            # We implement as:  g_linearized == sl_eq
            # Each element of noncvx_eq_lin is already "S[k+1] == rhs_lin"
            # so sl_eq captures the violation from the perspective of penalty
            # We instead formulate: impose the linear equality, track violation
            # via the penalty on the *nonlinear* residual g_k.
            # (The AL penalty acts on g_k evaluated at the solution, not on a slack.)
            # → Just add the linearized equalities as hard constraints.
            constraints += noncvx_eq_lin

        # Linearized nonconvex inequalities
        noncvx_ineq_lin = pr.noncvx_ineq_relaxed(cvx_vars, ref)
        constraints += noncvx_ineq_lin

        # Convexified inexact inequalities (obstacles)
        inexact_ineq = pr.convexified_inexact_ineq(cvx_vars, ref)
        constraints += inexact_ineq

        # Trust region
        constraints += self._trust_region_constraints(cvx_vars, ref, it.r)

        # ---- Objective  =  J0  +  AL penalty (acts on *predicted* slack) ----
        J0_expr = pr.objective(cvx_vars)

        # AL penalty on nonconvex equality residuals (predicted at linear model)
        penalty_expr = cp.Constant(0)
        # We use penalty on the actual residuals post-solve; the subproblem
        # minimises J0 subject to linearized constraints (exact AL update outside).

        objective = cp.Minimize(J0_expr)
        prob_cvx  = cp.Problem(objective, constraints)

        # ---- Solve ----
        solver_kwargs: dict[str, Any] = {
            "solver": p.solver,
            "verbose": p.solver_verbose,
        }
        if p.solver == "CLARABEL":
            solver_kwargs["eps_abs"] = p.solver_eps_abs
            solver_kwargs["eps_rel"] = p.solver_eps_rel
        elif p.solver == "MOSEK":
            solver_kwargs["mosek_params"] = {
                "MSK_DPAR_INTPNT_CO_TOL_PFEAS": p.solver_eps_abs,
                "MSK_DPAR_INTPNT_CO_TOL_DFEAS": p.solver_eps_abs,
            }

        try:
            prob_cvx.solve(**solver_kwargs)
        except cp.SolverError as e:
            return str(e), None, np.inf, np.array([]), np.array([]), np.array([])

        status = prob_cvx.status
        if prob_cvx.value is None or status in ("infeasible", "unbounded"):
            return status, None, np.inf, np.array([]), np.array([]), np.array([])

        J0_val = float(J0_expr.value)
        vals   = pr.extract_vals(cvx_vars)

        # Evaluate actual nonconvex residuals at the new solution
        g_vals = pr.noncvx_eq(vals)
        h_vals = pr.noncvx_ineq(vals)

        return (
            status,
            vals,
            J0_val,
            g_vals,                          # "slack_eq"    (actual residuals)
            h_vals,                          # "slack_ineq"  (actual residuals)
            np.array([]),                    # slack_inexact (not separately tracked)
        )

    # ================================================================== #
    # Trust region
    # ================================================================== #

    def _trust_region_constraints(
        self,
        cvx_vars: dict,
        ref: dict[str, np.ndarray],
        r: float,
    ) -> list[cp.Constraint]:
        """
        Physics-informed trust region for SqrtQRCovarianceSteering:
          For each k:  ||A[k](S[k] - S_ref[k]) + B[k](L[k] - L_ref[k])||_F <= r

        Falls back to standard variable-space trust region for other problems.
        """
        pr = self.prob
        p  = self.params

        if isinstance(pr, SqrtQRCovarianceSteering):
            return self._physics_tr(cvx_vars, ref, r, pr)
        else:
            return self._standard_tr(cvx_vars, ref, r)

    def _physics_tr(
        self,
        cvx_vars: dict,
        ref: dict[str, np.ndarray],
        r: float,
        pr: SqrtQRCovarianceSteering,
    ) -> list[cp.Constraint]:
        """Physics-informed trust region on linearised covariance propagation."""
        S_var, L_var = cvx_vars["S"], cvx_vars["L"]
        S_ref, L_ref = ref["S"], ref["L"]
        constraints = []

        for k in range(pr.N):
            dS = S_var[k] - S_ref[k]
            dL = L_var[k] - L_ref[k]
            # Scaled perturbation in propagation space
            delta_prop = pr.A[k] @ dS + pr.B[k] @ dL    # (nx, nx) CVXPY
            if p.trust_region_norm == np.inf:
                constraints.append(cp.norm_inf(cp.vec(delta_prop)) <= r)
            else:
                constraints.append(cp.norm(delta_prop, "fro") <= r)

        return constraints

    def _standard_tr(
        self,
        cvx_vars: dict,
        ref: dict[str, np.ndarray],
        r: float,
    ) -> list[cp.Constraint]:
        """Standard Euclidean trust region on all variables."""
        p = self.params
        constraints = []
        for key in cvx_vars:
            if key not in ref:
                continue
            var_list = cvx_vars[key]
            ref_arr  = ref[key]
            for k, var in enumerate(var_list):
                delta = var - ref_arr[k]
                if p.trust_region_norm == np.inf:
                    constraints.append(cp.norm_inf(cp.vec(delta)) <= r)
                else:
                    constraints.append(cp.norm(delta, "fro") <= r)
        return constraints

    # ================================================================== #
    # Penalty function (Augmented Lagrangian)
    # ================================================================== #

    def _penalty(
        self,
        g: np.ndarray,
        h: np.ndarray,
        slack_eq: np.ndarray,
        slack_ineq: np.ndarray,
        it: SCPIter,
    ) -> float:
        """
        AL penalty:
          P = lam' g + (w/2)||g||^2 + mu' max(h,0) + (w/2)||max(h,0)||^2
        """
        P = 0.0
        if len(g) > 0:
            P += float(it.lam @ g) + (it.w / 2) * float(g @ g)
        if len(h) > 0:
            h_pos = np.maximum(h, 0.0)
            P += float(it.mu_al @ h_pos) + (it.w / 2) * float(h_pos @ h_pos)
        return P

    # ================================================================== #
    # Helpers
    # ================================================================== #

    def _eval_objective_at_ref(self, ref: dict[str, np.ndarray]) -> float:
        """Evaluate the objective numerically at the reference point."""
        pr = self.prob
        if not isinstance(pr, SqrtQRCovarianceSteering):
            return 0.0

        S, L = ref["S"], ref["L"]
        J = 0.0
        for k in range(pr.N):
            J += np.trace(pr.Q @ S[k] @ S[k].T) + np.trace(pr.R @ L[k] @ L[k].T)
            if pr.use_mean and "mu" in ref:
                J += float(ref["mu"][k] @ pr.Q @ ref["mu"][k])
                J += float(ref["v"][k]  @ pr.R @ ref["v"][k])
        J += np.trace(pr.Q @ S[pr.N] @ S[pr.N].T)
        if pr.use_mean and "mu" in ref:
            J += float(ref["mu"][pr.N] @ pr.Q @ ref["mu"][pr.N])
        return J

    def _variable_change_norm(
        self, vals: dict[str, np.ndarray], ref: dict[str, np.ndarray]
    ) -> float:
        """||z_new - z_old||_inf across all variables."""
        diffs = []
        for key in vals:
            if key in ref:
                diffs.append(np.abs(vals[key] - ref[key]).ravel())
        if not diffs:
            return 0.0
        return float(np.max(np.concatenate(diffs)))

    def _warm_start(self, cvx_vars: dict, ref: dict[str, np.ndarray]) -> None:
        """Set CVXPY warm-start values from reference."""
        for key, var_list in cvx_vars.items():
            if key not in ref:
                continue
            ref_arr = ref[key]
            for k, var in enumerate(var_list):
                try:
                    var.value = ref_arr[k]
                except Exception:
                    pass

    # ================================================================== #
    # Printing
    # ================================================================== #

    def _print_header(self) -> None:
        print(
            f"{'iter':>5} | {'J0':>12} | {'chi':>10} | {'dJ':>10} | "
            f"{'rho':>7} | {'r':>8} | {'w':>8} | ref | mul"
        )
        print("-" * 82)

    def _print_iter(
        self, k, J0, chi, dJ, rho, r, w, ref_upd, mul_upd
    ) -> None:
        print(
            f"{k:5d} | {J0:12.4e} | {chi:10.3e} | {dJ:10.3e} | "
            f"{rho:7.3f} | {r:8.2e} | {w:8.2e} | "
            f"{'Y' if ref_upd else 'N':^3} | {'Y' if mul_upd else 'N':^3}"
        )
