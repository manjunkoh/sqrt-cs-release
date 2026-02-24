"""
Basic unit tests for the sqrt QR covariance steering components.

Run:  pytest tests/ -v
"""

import numpy as np
import pytest
from sqrt_cs.utils.qr_utils import economy_qr_pos_diag, d_QR, cholesky_pos_diag
from sqrt_cs.systems.double_integrator import build_double_integrator, default_boundary_conditions
from sqrt_cs.problems.sqrt_qr_cs import SqrtQRCovarianceSteering


class TestQRUtils:
    def test_positive_diagonal(self):
        rng = np.random.default_rng(0)
        M = rng.standard_normal((8, 4))
        Q, R = economy_qr_pos_diag(M)
        assert np.all(np.diag(R) >= 0)
        assert np.allclose(Q @ R, M, atol=1e-12)

    def test_orthonormality(self):
        rng = np.random.default_rng(1)
        M = rng.standard_normal((6, 4))
        Q, R = economy_qr_pos_diag(M)
        assert np.allclose(Q.T @ Q, np.eye(4), atol=1e-12)

    def test_d_QR_finite_difference(self):
        """Verify d_QR against numerical finite difference."""
        rng = np.random.default_rng(42)
        m, n = 8, 4
        M  = rng.standard_normal((m, n))
        dM = rng.standard_normal((m, n))
        eps = 1e-6

        Q0, R0 = economy_qr_pos_diag(M)
        _, R1   = economy_qr_pos_diag(M + eps * dM)

        dR_fd  = (R1 - R0) / eps
        _, dR  = d_QR(M, dM, Q=Q0, R=R0)

        assert np.allclose(dR, dR_fd, atol=1e-5), \
            f"Max error: {np.abs(dR - dR_fd).max():.2e}"


class TestProblemSetup:
    def setup_method(self):
        A, B, G = build_double_integrator(dt=1.0, d=2)   # smaller for speed
        P0, Pf, _, _ = default_boundary_conditions(d=2)
        self.prob = SqrtQRCovarianceSteering(
            N=5, A=A, B=B, G=G, P0=P0, Pf=Pf,
            Q=np.eye(4), R=np.eye(2),
        )

    def test_initial_ref(self):
        ref = self.prob.initial_ref()
        assert "S" in ref and "L" in ref
        assert ref["S"].shape == (6, 4, 4)
        assert ref["L"].shape == (5, 2, 4)

    def test_noncvx_eq_shape(self):
        ref = self.prob.initial_ref()
        g = self.prob.noncvx_eq(ref)
        # N * nx * nx residuals
        assert g.shape == (5 * 4 * 4,)

    def test_build_M_shape(self):
        S = np.eye(4) * 0.1
        L = np.zeros((2, 4))
        M = self.prob._build_M(S, L, 0)
        # (nx + nw, nx) = (4+2, 4) = (6, 4)
        assert M.shape == (6, 4)

    def test_covariance_propagation_loss_zero_at_qr_solution(self):
        """If S[k+1] is set to R.T from the QR, loss should be near zero."""
        from sqrt_cs.utils.qr_utils import economy_qr_pos_diag
        ref = self.prob.initial_ref()
        S, L = ref["S"], ref["L"]

        # Manually set S[k+1] = R[k].T for all k
        for k in range(self.prob.N):
            M = self.prob._build_M(S[k], L[k], k)
            _, R = economy_qr_pos_diag(M)
            S[k + 1] = R.T

        ref["S"] = S
        losses = self.prob.compute_covariance_propagation_loss(ref)
        assert np.all(losses < 1e-12), f"Losses: {losses}"
