"""
Unconstrained covariance steering — 3-D double integrator.

Reproduces the unconstrained benchmark from Section IV-A of the paper.
Run:
    python examples/example_unconstrained.py
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import time
import numpy as np
import matplotlib.pyplot as plt

from sqrt_cs.systems.double_integrator import build_double_integrator
from sqrt_cs.problems.sqrt_qr_cs import SqrtQRCovarianceSteering
from sqrt_cs.core.scp_params import SCPParams
from sqrt_cs.core.scvx_star import SCvxStar


def run(N: int = 20, verbose: bool = True):
    """
    Solve unconstrained CS for a 3-D double integrator with horizon N.

    Matches MATLAB horizon_size_scalability.m benchmark (Section IV-A of paper):
      - T_total = 3.0 s, dt = T_total / N
      - nw = nx = 6  (noise injected into all states)
      - G = sqrt(q * dt) * I_6  with q = 0.05
      - P0 = I_6,  Pf = 0.5 * I_6
      - Q = 0.1 * I_6,  R = I_3

    Parameters
    ----------
    N       : horizon length
    verbose : print convergence table

    Returns
    -------
    problem : solved SqrtQRCovarianceSteering instance
    solver  : SCvxStar instance (contains iteration history)
    """
    # ---- System (MATLAB-matching parameters) ----
    d = 3
    T_total = 3.0
    dt = T_total / N
    q = 0.05   # process noise spectral density

    A, B, _ = build_double_integrator(dt=dt, d=d)
    nx, nu = A.shape[0], B.shape[1]

    # Process noise: G = sqrt(q*dt) * I_{nx}  (nw = nx = 6, all-state noise)
    nw = nx
    G = np.sqrt(q * dt) * np.eye(nw)

    # Boundary conditions
    P0 = np.eye(nx)           # Sigma0 = I_6
    Pf = 0.5 * np.eye(nx)    # SigmaN = 0.5 * I_6

    # Cost matrices
    Q = 0.1 * np.eye(nx)
    R = np.eye(nu)

    # ---- Problem ----
    problem = SqrtQRCovarianceSteering(
        N=N,
        A=A, B=B, G=G,
        P0=P0, Pf=Pf,
        Q=Q, R=R,
        # No mean constraints → no mean variables (unconstrained case)
    )

    # ---- Solver params ----
    params = SCPParams(
        r_init=0.1,
        r_max=1.0,
        w_init=100.0,
        tol_opt=1e-5,
        tol_feas=1e-5,
        k_max=500,
        solver="CLARABEL",
        verbose=verbose,
    )

    # ---- Solve ----
    solver = SCvxStar(problem, params)
    t0 = time.perf_counter()
    sol = solver.solve()
    elapsed = time.perf_counter() - t0

    print(f"\nSolve time : {elapsed:.3f} s")
    print(f"Iterations : {len(solver.iter_state.hist_J0)}")

    if problem.P_sol is not None:
        print(f"Final cost : {solver.iter_state.hist_J0[-1]:.6e}")

        # ---- Covariance propagation loss ----
        loss = problem.compute_covariance_propagation_loss(sol)
        print(f"Max cov propagation loss : {loss.max():.2e}")
        print(f"Mean cov propagation loss: {loss.mean():.2e}")

    return problem, solver


def plot_results(problem: SqrtQRCovarianceSteering, solver: SCvxStar):
    """Plot convergence history and covariance trajectory."""
    it = solver.iter_state
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    # ---- Objective ----
    axes[0].semilogy(it.hist_J0, color="steelblue")
    axes[0].set_xlabel("SCP iteration")
    axes[0].set_ylabel("Objective J0")
    axes[0].set_title("Objective convergence")
    axes[0].grid(True, alpha=0.4)

    # ---- Feasibility ----
    axes[1].semilogy(it.hist_chi, color="darkorange")
    axes[1].set_xlabel("SCP iteration")
    axes[1].set_ylabel("Feasibility chi")
    axes[1].set_title("Constraint violation")
    axes[1].grid(True, alpha=0.4)

    # ---- State std devs over horizon ----
    if problem.P_sol is not None:
        N = problem.N
        sigmas = np.array([np.sqrt(np.diag(problem.P_sol[k])) for k in range(N + 1)])
        for i in range(problem.nx):
            axes[2].plot(sigmas[:, i], label=f"x{i+1}")
        axes[2].set_xlabel("Time step k")
        axes[2].set_ylabel("Standard deviation")
        axes[2].set_title("State std dev trajectory")
        axes[2].legend(fontsize=7, ncol=2)
        axes[2].grid(True, alpha=0.4)

    plt.tight_layout()
    plt.savefig("unconstrained_results.png", dpi=150)
    print("Saved unconstrained_results.png")
    plt.show()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=20, help="Horizon length")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()

    problem, solver = run(N=args.N)
    if not args.no_plot:
        plot_results(problem, solver)
