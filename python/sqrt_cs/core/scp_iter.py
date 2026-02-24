"""
SCPIter — mutable state carried across SCP iterations.

Stores Lagrange multipliers, penalty weights, trust region radius,
and per-iteration history for diagnostics / plotting.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np


@dataclass
class SCPIter:
    """
    Mutable state for one SCP iteration.

    Parameters
    ----------
    n_eq   : number of nonconvex equality residuals   (length of g)
    n_ineq : number of nonconvex inequality residuals (length of h)
    n_inexact_ineq : number of convexified-inexact inequality residuals
    r_init : initial trust-region radius
    w_init : initial penalty weight
    """
    n_eq:   int
    n_ineq: int
    n_inexact_ineq: int
    r_init: float = 0.1
    w_init: float = 100.0

    # Augmented Lagrangian multipliers  (initialised to zero)
    lam:    np.ndarray = field(init=False)   # shape (n_eq,)   equality
    mu_al:  np.ndarray = field(init=False)   # shape (n_ineq,) inequality (>=0)

    # Penalty weight and stationarity tolerance
    w:     float = field(init=False)
    delta: float = field(init=False)

    # Trust region radius
    r: float = field(init=False)

    # ------------------------------------------------------------------ #
    # Per-iteration history (appended by SCvxStar)
    # ------------------------------------------------------------------ #
    hist_J0:        list[float] = field(default_factory=list)
    hist_chi:       list[float] = field(default_factory=list)   # feasibility
    hist_delta_J:   list[float] = field(default_factory=list)   # actual improvement
    hist_delta_L:   list[float] = field(default_factory=list)   # predicted improvement
    hist_rho:       list[float] = field(default_factory=list)   # ratio
    hist_r:         list[float] = field(default_factory=list)   # TR radius
    hist_w:         list[float] = field(default_factory=list)   # penalty weight
    hist_ref_updated: list[bool] = field(default_factory=list)
    hist_mul_updated: list[bool] = field(default_factory=list)

    def __post_init__(self):
        self.lam   = np.zeros(self.n_eq)
        self.mu_al = np.zeros(self.n_ineq)
        self.w     = self.w_init
        self.delta = np.inf
        self.r     = self.r_init

    def record(
        self,
        J0: float,
        chi: float,
        delta_J: float,
        delta_L: float,
        rho: float,
        ref_updated: bool,
        mul_updated: bool,
    ) -> None:
        self.hist_J0.append(J0)
        self.hist_chi.append(chi)
        self.hist_delta_J.append(delta_J)
        self.hist_delta_L.append(delta_L)
        self.hist_rho.append(rho)
        self.hist_r.append(self.r)
        self.hist_w.append(self.w)
        self.hist_ref_updated.append(ref_updated)
        self.hist_mul_updated.append(mul_updated)
