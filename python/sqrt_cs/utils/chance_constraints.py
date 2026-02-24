"""
Chance constraint data structures.

A chance constraint is a probabilistic safety requirement of the form

    P( event ) >= 1 - p

where p is a small violation probability (e.g. 0.001).

Two types are supported:

  AffineCC   :  P( alpha' x  <=  beta )  >= 1 - p
  CircularObstacle : P( ||pos - center||_2 >= radius ) >= 1 - p

Both are reformulated into second-order cone constraints in terms of
the Cholesky factor S_k (state covariance sqrt) and mean mu_k.

Reference: Section III-B of the paper.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.stats import norm, chi2


@dataclass
class AffineCC:
    """
    Affine chance constraint:  P( alpha' x <= beta ) >= 1 - p

    SOCP reformulation (using S_k = chol(P_k)):
        alpha' mu_k + z_p * ||S_k' alpha||_2 <= beta
    where z_p = norm.ppf(1 - p).

    Attributes
    ----------
    alpha : (nx,)  direction vector
    beta  : float  upper bound
    p     : float  violation probability  (default 0.001 → 3-sigma-ish)
    nodes : list[int] | None
        Time indices at which to impose the constraint.
        None means all nodes 0..N (inclusive).
    """
    alpha: np.ndarray
    beta: float
    p: float = 0.001
    nodes: list[int] | None = None

    def __post_init__(self):
        self.alpha = np.asarray(self.alpha, dtype=float)
        self._z = float(norm.ppf(1.0 - self.p))

    @property
    def z(self) -> float:
        """Inverse normal CDF at 1-p."""
        return self._z


@dataclass
class NormCC:
    """
    Norm chance constraint:  P( ||x||_2 <= gamma ) >= 1 - p

    SOCP reformulation:
        ||mu_k||_2 + q_p * ||S_k||_F <= gamma
    where q_p = sqrt( chi2.ppf(1-p, df=nx) ).

    Attributes
    ----------
    gamma : float  radius bound
    nx    : int    state dimension (degrees of freedom for chi-squared)
    p     : float  violation probability
    nodes : list[int] | None
    """
    gamma: float
    nx: int
    p: float = 0.001
    nodes: list[int] | None = None

    def __post_init__(self):
        self._q = float(np.sqrt(chi2.ppf(1.0 - self.p, df=self.nx)))

    @property
    def q(self) -> float:
        """sqrt of chi-squared quantile at 1-p with nx dof."""
        return self._q


@dataclass
class CircularObstacle:
    """
    Circular (2-D) / spherical obstacle avoidance chance constraint:
        P( ||pos - center||_2 >= radius ) >= 1 - p

    The constraint is nonconvex but is linearized around a reference mean
    mu_ref at each SCP iteration into a SOCP half-space.

    Linearized form (Section III-C):
        a' mu_k + z_p * ||S_k[pos_idx, :] ' a||  >=  radius + a' center
    where a = (mu_ref[pos_idx] - center) / ||mu_ref[pos_idx] - center||

    Attributes
    ----------
    center   : (d,)  obstacle center (in position subspace)
    radius   : float keep-out radius
    pos_idx  : list[int]  indices of position states in full state vector
    p        : float  collision probability tolerance
    nodes    : list[int] | None
    """
    center: np.ndarray
    radius: float
    pos_idx: list[int]
    p: float = 0.001
    nodes: list[int] | None = None

    def __post_init__(self):
        self.center = np.asarray(self.center, dtype=float)
        self._z = float(norm.ppf(1.0 - self.p))

    @property
    def z(self) -> float:
        return self._z

    def linearize(self, mu_ref: np.ndarray) -> tuple[np.ndarray, float]:
        """
        Compute the linearization normal vector and offset at mu_ref.

        Parameters
        ----------
        mu_ref : (nx,) reference mean at a given time step

        Returns
        -------
        a      : (d,)  unit normal pointing away from obstacle center
        offset : float  RHS constant  (radius - a' center ... rearranged)
        """
        pos_ref = mu_ref[self.pos_idx]
        diff = pos_ref - self.center
        dist = np.linalg.norm(diff)
        if dist < 1e-10:
            # degenerate: perturb slightly
            diff = np.ones_like(self.center) / np.sqrt(len(self.center))
            dist = 1.0
        a = diff / dist
        return a, self.radius
