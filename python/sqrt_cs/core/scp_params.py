"""
SCP (Sequential Convex Programming) hyperparameters.

Mirrors SCPParams.m from the MATLAB codebase.
All defaults reproduce the paper's numerical results.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Literal


@dataclass
class SCPParams:
    # ------------------------------------------------------------------ #
    # Trust region
    # ------------------------------------------------------------------ #
    r_init: float = 0.1          # initial radius
    r_min:  float = 1e-6         # minimum radius
    r_max:  float = 1.0          # maximum radius
    trust_region_norm: float = np.inf  # norm order: np.inf or 2

    # ------------------------------------------------------------------ #
    # Augmented Lagrangian penalty
    # ------------------------------------------------------------------ #
    w_init:    float = 100.0     # initial penalty weight
    w_max:     float = 1e8       # maximum weight
    beta:      float = 2.0       # weight increase multiplier
    gamma_al:  float = 0.9       # stationarity tolerance multiplier
    w_inexact: float = 1000.0    # weight for convexified-inexact constraints

    penalty_method: Literal["AL", "L1", "ALwithL1"] = "AL"

    # ------------------------------------------------------------------ #
    # Convergence tolerances
    # ------------------------------------------------------------------ #
    tol_opt:    float = 1e-5    # |delta_J| <= tol_opt
    tol_feas:   float = 1e-5    # ||[g; max(h,0)]|| <= tol_feas
    tol_change: float = np.inf  # ||Delta z||_inf <= tol_change
    k_max:      int   = 2000    # maximum SCP iterations

    # ------------------------------------------------------------------ #
    # Step acceptance ratios (exact linearization)
    # ------------------------------------------------------------------ #
    rho0: float = 0.0    # accept if rho >= rho0
    rho1: float = 0.25   # shrink TR if rho < rho1
    rho2: float = 0.70   # expand TR if rho >= rho2

    # ------------------------------------------------------------------ #
    # Trust region update multipliers
    # ------------------------------------------------------------------ #
    alpha1: float = 2.0   # divisor  when shrinking
    alpha2: float = 3.0   # multiplier when expanding

    # ------------------------------------------------------------------ #
    # CVXPY solver
    # ------------------------------------------------------------------ #
    solver: str = "CLARABEL"     # "CLARABEL" (free) or "MOSEK" (license)
    solver_verbose: bool = False
    solver_eps_abs: float = 1e-8
    solver_eps_rel: float = 1e-8

    # ------------------------------------------------------------------ #
    # Misc
    # ------------------------------------------------------------------ #
    verbose: bool = True         # print convergence table
    linearization: Literal["exact", "inexact"] = "exact"
