"""
3-D Double Integrator  (matches paper's unconstrained benchmark, Section IV-A).

State:   x = [pos_x, pos_y, pos_z, vel_x, vel_y, vel_z]    (nx = 6)
Control: u = [acc_x, acc_y, acc_z]                          (nu = 3)
Noise:   w = [noise_x, noise_y, noise_z]                    (nw = 3)

Continuous:
    p_dot = v,   v_dot = u + w

Discrete (ZOH, timestep dt):
    A = [[I, dt*I],
         [0,    I]]
    B = [[0.5*dt^2 * I],
         [dt * I      ]]
    G = [[0.5*dt^2 * I],
         [dt * I      ]]
"""

from __future__ import annotations

import numpy as np


def build_double_integrator(
    dt: float = 1.0, d: int = 3
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build LTI matrices for a d-dimensional double integrator.

    Parameters
    ----------
    dt : timestep
    d  : spatial dimension  (default 3 → 3-D, nx=6, nu=3, nw=3)

    Returns
    -------
    A : (2d, 2d)
    B : (2d, d)
    G : (2d, d)
    """
    I  = np.eye(d)
    Z  = np.zeros((d, d))

    A = np.block([[I, dt * I],
                  [Z,      I]])

    B = np.block([[0.5 * dt**2 * I],
                  [dt * I         ]])

    G = B.copy()   # same noise shaping as control

    return A, B, G


def default_boundary_conditions(
    d: int = 3,
    sigma0: float = 0.1,
    sigmaf: float = 0.01,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Default initial / terminal covariances and means for benchmarking.

    P0 = sigma0^2 * I_{2d}
    Pf = sigmaf^2 * I_{2d}
    mu0 = zeros(2d)
    muf = zeros(2d)   (mean steering is unconstrained by default)
    """
    nx = 2 * d
    P0 = sigma0**2 * np.eye(nx)
    Pf = sigmaf**2 * np.eye(nx)
    mu0 = np.zeros(nx)
    muf = np.zeros(nx)
    return P0, Pf, mu0, muf
