"""
QR decomposition utilities for square root covariance steering.

Key functions:
- economy_qr_pos_diag : thin QR with positive diagonal on R (unique decomposition)
- d_QR                : differential of the QR map  d(Q,R) = dQR(M, dM)
                        used to linearize the covariance dynamics in the SCP loop

Reference:
  Z. Lin, "Riemannian Geometry of Symmetric Positive Definite Matrices via
  Cholesky Decomposition", SIAM J. Matrix Anal. Appl., 2019.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import solve_triangular


def economy_qr_pos_diag(M: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Thin (economy) QR decomposition M = Q R with the sign convention
    that diagonal entries of R are non-negative (making it unique when M
    has full column rank).

    Parameters
    ----------
    M : (m, n)  matrix with m >= n

    Returns
    -------
    Q : (m, n)  orthonormal columns
    R : (n, n)  upper triangular, diag >= 0
    """
    Q, R = np.linalg.qr(M, mode="reduced")
    signs = np.sign(np.diag(R))
    signs[signs == 0] = 1.0          # treat exact zeros as positive
    Q = Q * signs[np.newaxis, :]
    R = R * signs[:, np.newaxis]
    return Q, R


def d_QR(
    M: np.ndarray,
    dM: np.ndarray,
    Q: np.ndarray | None = None,
    R: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Differential of the economy QR decomposition with positive diagonal.

    Given  M = Q R  and a perturbation  dM, returns  (dQ, dR)  such that
        M + eps*dM  ≈  (Q + eps*dQ)(R + eps*dR)   as eps → 0.

    Algorithm (Lin 2019, Eq. A.1):
        Omega = Q' dM R^{-1} - R^{-T} dM' Q         (skew-symmetric)
        dR    = (triu(Q' dM) + tril(Omega, -1) R      -- compact form below
        dQ    = (dM - Q dR) R^{-1}

    The standard efficient form avoids forming R^{-1} explicitly by using
    triangular solves.

    Parameters
    ----------
    M  : (m, n)
    dM : (m, n)  perturbation direction
    Q  : (m, n)  pre-computed Q factor  (recomputed if None)
    R  : (n, n)  pre-computed R factor  (recomputed if None)

    Returns
    -------
    dQ : (m, n)
    dR : (n, n)
    """
    if Q is None or R is None:
        Q, R = economy_qr_pos_diag(M)

    n = R.shape[0]

    # C = Q' dM  (n x n)
    C = Q.T @ dM                                           # (n, n)

    # T = C R^{-1}  via triangular solve:  R^T T^T = C^T
    T = solve_triangular(R.T, C.T, lower=True).T          # T = C R^{-1}  (n, n)

    # Omega = tril(T, -1) - tril(T, -1)'   (skew-symmetric)
    T_lo  = np.tril(T, -1)
    Omega = T_lo - T_lo.T                                  # (n, n)

    # dR = C - Omega @ R   (upper triangular by construction)
    dR = C - Omega @ R                                     # (n, n)

    # dQ = (dM - Q dR) R^{-1}   via triangular solve
    dQ_R = dM - Q @ dR                                     # (m, n)
    dQ = solve_triangular(R.T, dQ_R.T, lower=True).T      # (m, n)

    return dQ, dR


def cholesky_pos_diag(P: np.ndarray) -> np.ndarray:
    """
    Lower Cholesky factor of P with positive diagonal.
    P = S S'   where S is lower triangular, diag(S) > 0.
    """
    S = np.linalg.cholesky(P)
    # np.linalg.cholesky already returns lower triangular with positive diagonal
    return S


def reconstruct_covariance(S: np.ndarray) -> np.ndarray:
    """P = S S'  (S lower triangular)."""
    return S @ S.T
