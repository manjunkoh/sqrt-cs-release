"""
Square Root QR Covariance Steering  (linear time-varying systems).

Implements the formulation from:
  "Square Root-Factorized Covariance Steering", 2025.

Decision variables
------------------
S[k]  : (nx, nx)  lower-triangular Cholesky factor of state covariance P[k]
L[k]  : (nu, nx)  feedback gain factor  (K[k] = L[k] @ inv(S[k]))
mu[k] : (nx,)     state mean  (optional; required when mean boundary conditions given)
v[k]  : (nu,)     control mean (optional; required when mean boundary conditions given)

Covariance dynamics (nonconvex)
---------
Build  M[k] = vstack( A[k] S[k] + B[k] L[k],  G[k] )        shape (nx+nw, nx)
Then   QR factorize:  M[k] = Q[k] R[k]
Set    S[k+1] = R[k].T                                         (lower triangular)

This is the key nonlinear step that the SCP linearizes.

Extension hook
--------------
Subclass and override  `_build_M`  +  `_noncvx_residual_k`  to plug in
sigma-point (unscented) propagation for nonlinear systems (Phase 2, Option A).
"""

from __future__ import annotations

from typing import Any

import cvxpy as cp
import numpy as np

from sqrt_cs.core.scp_problem import SCPProblem
from sqrt_cs.utils.qr_utils import economy_qr_pos_diag, d_QR, cholesky_pos_diag
from sqrt_cs.utils.chance_constraints import AffineCC, NormCC, CircularObstacle


class SqrtQRCovarianceSteering(SCPProblem):
    """
    Square root QR covariance steering for LTV systems.

    Parameters
    ----------
    N      : horizon length
    A, B, G: system matrices, each shape (N, nx, nx/nu/nw) or (nx, nx/nu/nw) if LTI
    P0, Pf : initial / terminal covariance  (nx, nx)
    Q, R   : state / control cost matrices  (nx, nx) / (nu, nu)
    mu0, muf : initial / terminal mean  (nx,)  — if None, mean vars are not created
    chance_constraints_state   : list of AffineCC / NormCC on state
    chance_constraints_control : list of AffineCC on control
    obstacles : list of CircularObstacle
    init_S : (N+1, nx, nx) initial guess for S  — defaults to chol(P0) broadcast
    init_L : (N, nu, nx)   initial guess for L  — defaults to zeros
    """

    def __init__(
        self,
        N: int,
        A: np.ndarray,
        B: np.ndarray,
        G: np.ndarray,
        P0: np.ndarray,
        Pf: np.ndarray,
        Q: np.ndarray | None = None,
        R: np.ndarray | None = None,
        mu0: np.ndarray | None = None,
        muf: np.ndarray | None = None,
        chance_constraints_state: list | None = None,
        chance_constraints_control: list | None = None,
        obstacles: list[CircularObstacle] | None = None,
        init_S: np.ndarray | None = None,
        init_L: np.ndarray | None = None,
        init_mu: np.ndarray | None = None,
        init_v:  np.ndarray | None = None,
    ):
        self.N  = N
        self.nx = P0.shape[0]
        self.nu = B.shape[-1]
        self.nw = G.shape[-1]

        # Broadcast LTI to LTV
        self.A = self._broadcast(A, (self.nx, self.nx), N)   # (N, nx, nx)
        self.B = self._broadcast(B, (self.nx, self.nu), N)   # (N, nx, nu)
        self.G = self._broadcast(G, (self.nx, self.nw), N)   # (N, nx, nw)

        self.P0 = np.asarray(P0, dtype=float)
        self.Pf = np.asarray(Pf, dtype=float)
        self.S0 = cholesky_pos_diag(self.P0)
        self.Sf = cholesky_pos_diag(self.Pf)

        nx, nu = self.nx, self.nu
        self.Q = np.eye(nx) if Q is None else np.asarray(Q, dtype=float)
        self.R = np.eye(nu) if R is None else np.asarray(R, dtype=float)

        # Mean tracking
        self.use_mean = (mu0 is not None) or (muf is not None)
        self.mu0 = np.zeros(nx) if mu0 is None else np.asarray(mu0, dtype=float)
        self.muf = np.zeros(nx) if muf is None else np.asarray(muf, dtype=float)

        # Constraints
        self.cc_state   = chance_constraints_state   or []
        self.cc_control = chance_constraints_control or []
        self.obstacles  = obstacles or []

        # Initial guesses
        if init_S is not None:
            self.init_S = init_S
        else:
            # Linearly interpolate from S0 to Sf to satisfy the terminal constraint
            S_interp = np.zeros((N + 1, nx, nx))
            for k in range(N + 1):
                alpha = k / N
                S_interp[k] = (1 - alpha) * self.S0 + alpha * self.Sf
            self.init_S = S_interp
        self.init_L  = init_L  if init_L  is not None else np.zeros((N, nu, nx))
        self.init_mu = init_mu if init_mu is not None else np.tile(self.mu0, (N + 1, 1))
        self.init_v  = init_v  if init_v  is not None else np.zeros((N, nu))

        # Populated after solve
        self.S_sol: np.ndarray | None = None   # (N+1, nx, nx)
        self.L_sol: np.ndarray | None = None   # (N, nu, nx)
        self.K_sol: np.ndarray | None = None   # (N, nu, nx) feedback gains
        self.P_sol: np.ndarray | None = None   # (N+1, nx, nx)
        self.mu_sol: np.ndarray | None = None  # (N+1, nx)
        self.v_sol:  np.ndarray | None = None  # (N, nu)

    # ================================================================== #
    # SCPProblem interface
    # ================================================================== #

    def define_vars(self) -> dict[str, Any]:
        N, nx, nu = self.N, self.nx, self.nu
        # S[k] is lower triangular — we use a full (nx, nx) variable and
        # enforce lower-triangular + positive diagonal in convex_ineq.
        vars = {
            "S":  [cp.Variable((nx, nx), name=f"S_{k}") for k in range(N + 1)],
            "L":  [cp.Variable((nu, nx), name=f"L_{k}") for k in range(N)],
        }
        if self.use_mean:
            vars["mu"] = [cp.Variable(nx, name=f"mu_{k}") for k in range(N + 1)]
            vars["v"]  = [cp.Variable(nu, name=f"v_{k}")  for k in range(N)]
        return vars

    def objective(self, vars: dict) -> cp.Expression:
        """
        LQG cost (matches MATLAB 'LQG'/'LQR'/'LQ' objective):
          J = sum_{k=0}^{N-1} [ trace(Q P[k]) + trace(R P_u[k]) ]
        where  P[k] = S[k] S[k]'  and  P_u[k] = L[k] L[k]' + (mean contribution).

        No terminal cost on S[N] — matches MATLAB SqrtQRCovarianceSteering.

        In CVXPY we use the identity  trace(Q S S') = ||chol(Q) S||_F^2.
        """
        S, L = vars["S"], vars["L"]
        N = self.N

        cost = cp.Constant(0)

        # Pre-compute Cholesky of cost matrices once
        Lq = np.linalg.cholesky(self.Q)   # Q = Lq Lq'
        Lr = np.linalg.cholesky(self.R)   # R = Lr Lr'

        for k in range(N):
            cost = cost + cp.sum_squares(Lq @ S[k]) + cp.sum_squares(Lr @ L[k])

        # Mean cost (if tracking means)
        if self.use_mean:
            mu, v = vars["mu"], vars["v"]
            for k in range(N):
                cost = cost + cp.quad_form(mu[k], self.Q) + cp.quad_form(v[k], self.R)
            cost = cost + cp.quad_form(mu[N], self.Q)

        return cost

    def convex_eq(self, vars: dict) -> list[cp.Constraint]:
        """
        Always-convex equality constraints:
          - Initial covariance:  S[0] = S0
          - Lower-triangular structure: upper triangle of S[k] = 0
          - Mean boundary conditions (if use_mean)
        """
        S = vars["S"]
        constraints = []

        # Initial covariance factor
        constraints.append(S[0] == self.S0)

        # Lower-triangular structure for all k
        nx = self.nx
        for k in range(self.N + 1):
            for i in range(nx):
                for j in range(i + 1, nx):
                    constraints.append(S[k][i, j] == 0)

        # Mean dynamics (linear — always convex)
        if self.use_mean:
            mu, v = vars["mu"], vars["v"]
            # Initial mean
            constraints.append(mu[0] == self.mu0)
            # Mean propagation:  mu[k+1] = A[k] mu[k] + B[k] v[k]
            for k in range(self.N):
                constraints.append(
                    mu[k + 1] == self.A[k] @ mu[k] + self.B[k] @ v[k]
                )

        return constraints

    def convex_ineq(self, vars: dict) -> list[cp.Constraint]:
        """
        Always-convex inequality constraints:
          - Positive diagonal on S[k]  (ensures positive definiteness)
          - Terminal covariance: ||chol(Pf)^{-1} S[N]||_2 <= 1  (PSD upper bound P[N] <= Pf)
            Matches MATLAB: norm(chol(P_f,'lower') \\ S[:,N], 2) <= 1
          - Chance constraints on state and control (convex in S, mu)
          - Terminal mean (if use_mean)
        """
        S, L = vars["S"], vars["L"]
        N, nx = self.N, self.nx
        constraints = []

        # Positive diagonal entries of S[k]
        for k in range(N + 1):
            for i in range(nx):
                constraints.append(S[k][i, i] >= 1e-8)

        # Terminal covariance: spectral norm inequality ||Sf^{-1} S[N]||_2 <= 1
        # Equivalent to P[N] = S[N]S[N]' <= Pf  (PSD inequality)
        # Sf = chol(Pf, 'lower') is stored as self.Sf
        Sf_inv = np.linalg.inv(self.Sf)
        constraints.append(cp.norm(Sf_inv @ S[N], 2) <= 1)

        # Terminal mean
        if self.use_mean:
            constraints.append(vars["mu"][N] == self.muf)

        # ---- State chance constraints ----
        for cc in self.cc_state:
            nodes = cc.nodes if cc.nodes is not None else list(range(N + 1))
            for k in nodes:
                constraints += self._affine_cc_constraint(S[k], vars, k, cc, "state")

        # ---- Control chance constraints ----
        for cc in self.cc_control:
            nodes = cc.nodes if cc.nodes is not None else list(range(N))
            for k in nodes:
                # P_u[k] = L[k] L[k]'  → chol factor of P_u is L[k]
                constraints += self._affine_cc_constraint(L[k], vars, k, cc, "control")

        return constraints

    def noncvx_eq(self, vals: dict[str, np.ndarray]) -> np.ndarray:
        """
        Evaluate QR covariance dynamics residuals numerically.

        For each k:  residual[k] = S[k+1] - R[k].T
        where  [A[k] S[k] + B[k] L[k]; G[k]] = Q[k] R[k]  (QR with pos diag).

        Returns a 1-D array of length  N * nx * nx.
        """
        S_vals = vals["S"]   # list of (nx, nx) arrays
        L_vals = vals["L"]   # list of (nu, nx) arrays
        residuals = []

        for k in range(self.N):
            M = self._build_M(S_vals[k], L_vals[k], k)
            _, R = economy_qr_pos_diag(M)
            residuals.append((S_vals[k + 1] - R.T).ravel())

        return np.concatenate(residuals)

    def noncvx_eq_relaxed(
        self, vars: dict, ref: dict[str, np.ndarray]
    ) -> list[cp.Expression]:
        """
        Linearized QR dynamics residuals:  S[k+1] - (R_ref[k].T + dR[k].T)

        where  dR[k]  is the differential of the QR map evaluated at
        (M_ref[k], dM[k])  with  dM[k] = (A[k] dS[k] + B[k] dL[k]).

        Returns a list of (nx, nx) CVXPY expressions, one per time step.
        Each expression equals zero when the linearized constraint is satisfied.
        """
        S_var, L_var = vars["S"], vars["L"]
        S_ref = ref["S"]    # (N+1, nx, nx)
        L_ref = ref["L"]    # (N, nu, nx)

        exprs = []
        for k in range(self.N):
            M_ref = self._build_M(S_ref[k], L_ref[k], k)
            Q_ref, R_ref = economy_qr_pos_diag(M_ref)

            # Differential of M w.r.t. (S[k], L[k]) perturbations
            dS = S_var[k] - S_ref[k]   # CVXPY expression
            dL = L_var[k] - L_ref[k]

            dM = self._build_dM_cvxpy(dS, dL, k)   # (nx+nw, nx) CVXPY expression

            # dR via d_QR — we need numpy version for the linear map
            # The differential is linear in dM, so we can apply it element-wise.
            dR_expr = self._dR_cvxpy(dM, Q_ref, R_ref)

            S_next_lin = R_ref.T + cp.reshape(dR_expr, (self.nx, self.nx), order='F').T

            # Return residual expression: S[k+1] - linearized_prediction
            residual = S_var[k + 1] - S_next_lin
            exprs.append(residual)

        return exprs

    def convexified_inexact_ineq(
        self, vars: dict, ref: dict[str, np.ndarray]
    ) -> list[cp.Constraint]:
        """
        Linearized obstacle-avoidance constraints (updated every SCP iteration).
        """
        if not self.obstacles or not self.use_mean:
            return []

        S_var = vars["S"]
        mu_var = vars["mu"]
        mu_ref = ref.get("mu")   # (N+1, nx)
        if mu_ref is None:
            return []

        constraints = []
        for obs in self.obstacles:
            nodes = obs.nodes if obs.nodes is not None else list(range(1, self.N + 1))
            for k in nodes:
                a, radius = obs.linearize(mu_ref[k])
                # Build selection matrix for position indices
                E = np.zeros((len(obs.pos_idx), self.nx))
                for i, idx in enumerate(obs.pos_idx):
                    E[i, idx] = 1.0

                # Linearized obstacle avoidance chance constraint (Sec III-C):
                #   P(a'x_pos >= a'center + radius) >= 1-p
                #   ↔  a'mu - z * ||S' a_full||_2 >= a'center + radius
                #   DCP form: z * norm(S' a_full) <= a'mu - a'center - radius
                # (norm <= affine is DCP; opposite of state upper-bound constraints)
                Ea = E.T @ a              # (nx,) unit normal in full-state space
                center_offset = float(a @ obs.center)  # a' center (constant)
                lhs_mean = Ea @ mu_var[k]
                lhs_cov  = cp.norm(S_var[k].T @ Ea, 2) * obs.z
                constraints.append(lhs_cov <= lhs_mean - obs.radius - center_offset)

        return constraints

    def postprocess(self, vals: dict[str, np.ndarray], ref: dict) -> None:
        """Extract P, K from S, L after convergence."""
        N, nx, nu = self.N, self.nx, self.nu
        S = vals["S"]   # (N+1, nx, nx)
        L = vals["L"]   # (N, nu, nx)

        self.S_sol = S
        self.L_sol = L
        self.P_sol = np.array([S[k] @ S[k].T for k in range(N + 1)])

        # K[k] = L[k] @ inv(S[k])  — triangular solve
        from scipy.linalg import solve_triangular
        self.K_sol = np.array([
            solve_triangular(S[k].T, L[k].T, lower=False).T
            for k in range(N)
        ])

        if self.use_mean:
            self.mu_sol = vals.get("mu")
            self.v_sol  = vals.get("v")

    # ================================================================== #
    # Reference packing / unpacking helpers (used by SCvxStar)
    # ================================================================== #

    def pack_ref(
        self,
        S: np.ndarray,
        L: np.ndarray,
        mu: np.ndarray | None = None,
        v:  np.ndarray | None = None,
    ) -> dict[str, np.ndarray]:
        ref = {"S": S.copy(), "L": L.copy()}
        if self.use_mean and mu is not None:
            ref["mu"] = mu.copy()
            ref["v"]  = v.copy() if v is not None else np.zeros((self.N, self.nu))
        return ref

    def initial_ref(self) -> dict[str, np.ndarray]:
        ref = {"S": self.init_S.copy(), "L": self.init_L.copy()}
        if self.use_mean:
            ref["mu"] = self.init_mu.copy()
            ref["v"]  = self.init_v.copy()
        return ref

    def extract_vals(self, vars: dict) -> dict[str, np.ndarray]:
        """Read CVXPY variable values into numpy arrays."""
        vals: dict[str, np.ndarray] = {}
        vals["S"] = np.array([vars["S"][k].value for k in range(self.N + 1)])
        vals["L"] = np.array([vars["L"][k].value for k in range(self.N)])
        if self.use_mean:
            vals["mu"] = np.array([vars["mu"][k].value for k in range(self.N + 1)])
            vals["v"]  = np.array([vars["v"][k].value  for k in range(self.N)])
        return vals

    def compute_covariance_propagation_loss(self, vals: dict[str, np.ndarray]) -> np.ndarray:
        """
        Per-step covariance propagation loss  (eq. 39 in paper):
          loss[k] = ||phi(K[k], P[k]) - P[k+1]||_F / ||P[k+1]||_F

        where phi is the exact covariance propagation map.
        """
        S = vals["S"]
        L = vals["L"]
        losses = []
        for k in range(self.N):
            M  = self._build_M(S[k], L[k], k)
            _, R = economy_qr_pos_diag(M)
            P_next_true  = R.T @ R
            P_next_stored = S[k + 1] @ S[k + 1].T
            num  = np.linalg.norm(P_next_true - P_next_stored, "fro")
            denom = np.linalg.norm(P_next_stored, "fro")
            losses.append(num / (denom + 1e-30))
        return np.array(losses)

    # ================================================================== #
    # Extension hook — override for nonlinear systems (Phase 2)
    # ================================================================== #

    def _build_M(self, S_k: np.ndarray, L_k: np.ndarray, k: int) -> np.ndarray:
        """
        Build the stacked matrix whose QR factorization gives S[k+1].

        Linear case:
            P_{k+1} = (A S + B L)(A S + B L)' + G G'
                    = M' M    where M = vstack([(A S + B L)', G'])
            shape of M: (nx + nw, nx)

        S_{k+1} = R'  from thin QR  M = Q R.

        Override this method (and _build_dM_cvxpy / _dR_cvxpy) to implement
        sigma-point (unscented) propagation for nonlinear systems.
        """
        top    = (self.A[k] @ S_k + self.B[k] @ L_k).T   # (nx, nx)
        bottom = self.G[k].T                               # (nw, nx)
        return np.vstack([top, bottom])                    # (nx+nw, nx)

    def _build_dM_cvxpy(self, dS: cp.Expression, dL: cp.Expression, k: int) -> cp.Expression:
        """
        CVXPY expression for the perturbation of M[k] w.r.t. (dS, dL).

        Linear case:  dM = [ (A[k] dS + B[k] dL)' ]    (nx, nx)
                           [           0           ]    (nw, nx)
        shape of dM: (nx + nw, nx)
        """
        nx, nw = self.nx, self.nw
        top    = (self.A[k] @ dS + self.B[k] @ dL).T    # (nx, nx) CVXPY
        bottom = cp.Constant(np.zeros((nw, nx)))
        return cp.vstack([top, bottom])                   # (nx+nw, nx)

    def _dR_cvxpy(
        self,
        dM: cp.Expression,
        Q_ref: np.ndarray,
        R_ref: np.ndarray,
    ) -> cp.Expression:
        """
        Linear map  dM → dR  via the d_QR differential.

        Since d_QR is a linear function of dM, we can apply it column-by-column
        and reconstruct as a CVXPY affine expression.

        Returns a CVXPY expression of shape (nx*nx,) representing dR.ravel().
        """
        nx = self.nx
        m, n = Q_ref.shape[0], R_ref.shape[0]

        # Represent dM as a CVXPY (m, n) expression.
        # We need to compute  dR  which is linear in dM.
        # Strategy: use numpy to compute the linear operator  T  such that
        #   dR.ravel() = T @ dM.ravel()
        # then express as  T @ cp.vec(dM).

        T = self._dR_linear_operator(Q_ref, R_ref)       # (n*n, m*n) numpy, F-order
        dM_vec = cp.vec(dM, order='F')                  # (m*n,) CVXPY, F-order
        return T @ dM_vec                               # (n*n,) CVXPY

    def _dR_linear_operator(self, Q: np.ndarray, R: np.ndarray) -> np.ndarray:
        """
        Build the  (n*n) × (m*n)  matrix representing the linear map  dM → dR.ravel().

        We use the standard basis:  apply d_QR to each standard basis vector
        in the dM space and collect the outputs.
        """
        m, n = Q.shape[0], R.shape[0]
        T = np.zeros((n * n, m * n))
        for col in range(m * n):
            e = np.zeros(m * n)
            e[col] = 1.0
            dM_basis = e.reshape(m, n, order='F')       # F-order: matches cp.vec(dM, order='F')
            _, dR = d_QR(None, dM_basis, Q=Q, R=R)
            T[:, col] = dR.ravel(order='F')             # F-order: matches cp.reshape(..., order='F')
        return T

    # ================================================================== #
    # Private helpers
    # ================================================================== #

    @staticmethod
    def _broadcast(
        mat: np.ndarray, shape_2d: tuple[int, int], N: int
    ) -> np.ndarray:
        """Broadcast a 2-D matrix to shape (N, *shape_2d) if needed."""
        mat = np.asarray(mat, dtype=float)
        if mat.ndim == 2:
            return np.tile(mat[np.newaxis], (N, 1, 1))
        assert mat.shape == (N, *shape_2d), f"Expected ({N}, {shape_2d[0]}, {shape_2d[1]}), got {mat.shape}"
        return mat

    def _affine_cc_constraint(
        self,
        chol_var: cp.Variable,
        vars: dict,
        k: int,
        cc: AffineCC | NormCC,
        kind: str,
    ) -> list[cp.Constraint]:
        """
        Convert an AffineCC / NormCC to CVXPY SOCP constraints.

        kind = 'state'   → chol_var is S[k],  mean is mu[k]
        kind = 'control' → chol_var is L[k],  mean is v[k]
        """
        constraints = []

        if isinstance(cc, AffineCC):
            a = cc.alpha
            if self.use_mean:
                mu_k = vars["mu"][k] if kind == "state" else vars["v"][k]
                mean_term = a @ mu_k
            else:
                mean_term = 0.0

            # ||chol_var' alpha||_2 * z + alpha' mu <= beta
            cone_arg = chol_var.T @ a
            constraints.append(
                mean_term + cc.z * cp.norm(cone_arg, 2) <= cc.beta
            )

        elif isinstance(cc, NormCC):
            if self.use_mean:
                mu_k = vars["mu"][k] if kind == "state" else vars["v"][k]
                constraints.append(
                    cp.norm(mu_k, 2) + cc.q * cp.normNuc(chol_var) <= cc.gamma
                    # normNuc of lower-triangular = Frobenius of S ≈ trace(sqrt(P))
                    # For tighter bound use Frobenius: ||S||_F >= ||mu||
                )
            else:
                constraints.append(cp.norm(chol_var, "fro") * cc.q <= cc.gamma)

        return constraints
