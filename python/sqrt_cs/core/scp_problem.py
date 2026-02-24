"""
Abstract base class for SCP sub-problems.

Every concrete covariance steering formulation (linear QR, nonlinear sigma-point, etc.)
inherits from SCPProblem and implements the methods below.

Design principle
----------------
The SCP solver (SCvxStar) only interacts with the problem through this interface.
This makes it trivial to swap in a new dynamics model (e.g. sigma-point propagation
for the nonlinear extension) without touching the solver.

Variable convention
-------------------
All CVXPY decision variables are stored in a dict  vars: dict[str, cp.Variable | list].
The reference (linearization point) is stored as    ref: dict[str, np.ndarray | list].
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

import cvxpy as cp
import numpy as np


class SCPProblem(ABC):
    """
    Abstract interface between SCvxStar and a specific trajectory-optimization problem.

    Subclasses must implement
    --------------------------
    define_vars()             → dict of CVXPY variables
    objective(vars)           → scalar CVXPY expression
    convex_eq(vars)           → list of CVXPY constraints  (always convex, affine)
    convex_ineq(vars)         → list of CVXPY constraints  (always convex)
    noncvx_eq(vals)           → np.ndarray  of residuals  g(z) = 0
    noncvx_eq_relaxed(vars, ref) → list of CVXPY constraints  (linearized g ≈ 0)

    Optionally override
    -------------------
    noncvx_ineq(vals)                      (default: empty)
    noncvx_ineq_relaxed(vars, ref)         (default: empty)
    convexified_inexact_ineq(vars, ref)    (default: empty — used for obstacles)
    vars_to_vec(vars)                      → 1-D np.ndarray  for trust region
    ref_to_vec(ref)                        → 1-D np.ndarray  for trust region
    postprocess(vars, ref)                 post-solution cleanup
    """

    # ------------------------------------------------------------------ #
    # Abstract — must implement
    # ------------------------------------------------------------------ #

    @abstractmethod
    def define_vars(self) -> dict[str, Any]:
        """Return a dict of CVXPY Variables (or lists thereof)."""

    @abstractmethod
    def objective(self, vars: dict) -> cp.Expression:
        """Return the scalar CVXPY objective expression."""

    @abstractmethod
    def convex_eq(self, vars: dict) -> list[cp.Constraint]:
        """
        Constraints that are always convex (affine equalities).
        Returned list is fixed across SCP iterations.
        """

    @abstractmethod
    def convex_ineq(self, vars: dict) -> list[cp.Constraint]:
        """
        Constraints that are always convex (e.g. terminal covariance LMI,
        positive diagonal of S, chance constraints convex in S).
        Returned list is fixed across SCP iterations.
        """

    @abstractmethod
    def noncvx_eq(self, vals: dict[str, np.ndarray]) -> np.ndarray:
        """
        Evaluate the nonconvex equality residuals  g(z)  numerically.
        vals  maps variable names to numpy arrays (evaluated solution).
        Returns a 1-D array; the solver drives this to zero.
        """

    @abstractmethod
    def noncvx_eq_relaxed(
        self, vars: dict, ref: dict[str, np.ndarray]
    ) -> list[cp.Constraint]:
        """
        First-order Taylor expansion of  g(z) ≈ 0  around  ref.
        Returns CVXPY constraints used inside the convex subproblem.
        """

    # ------------------------------------------------------------------ #
    # Optional — override if needed
    # ------------------------------------------------------------------ #

    def noncvx_ineq(self, vals: dict[str, np.ndarray]) -> np.ndarray:
        """Evaluate nonconvex inequality residuals  h(z) <= 0  numerically."""
        return np.array([])

    def noncvx_ineq_relaxed(
        self, vars: dict, ref: dict[str, np.ndarray]
    ) -> list[cp.Constraint]:
        """Linearized nonconvex inequality constraints."""
        return []

    def convexified_inexact_ineq(
        self, vars: dict, ref: dict[str, np.ndarray]
    ) -> list[cp.Constraint]:
        """
        Constraints that are nonconvex but admit a convex approximation
        (e.g. obstacle avoidance linearized as half-spaces).
        Updated every SCP iteration.
        """
        return []

    def vars_to_vec(self, vars: dict[str, np.ndarray]) -> np.ndarray:
        """
        Flatten evaluated variables to a single vector for trust-region distance.
        Default: concatenate all arrays in insertion order.
        Override for physics-informed trust regions.
        """
        pieces = []
        for v in vars.values():
            arr = np.asarray(v)
            pieces.append(arr.ravel())
        return np.concatenate(pieces) if pieces else np.array([])

    def postprocess(self, vars: dict, ref: dict) -> None:
        """Called by SCvxStar after convergence to extract the final solution."""

    def pre_iteration(self, k: int, ref: dict) -> None:
        """Hook called before each SCP iteration (e.g. to update parameters)."""

    def post_iteration(self, k: int, vars: dict, ref: dict) -> None:
        """Hook called after each SCP iteration."""
