"""
Constrained covariance steering — 2-D double integrator with obstacle avoidance.

Matches MATLAB horizon_size_scalability_single_obstacle.m (Section IV-B of paper):
  - 2-D double integrator: nx=4 [px,py,vx,vy], nu=2 [ax,ay]
  - Mean steering: [0,0,0,0] → [10,0,0,0]
  - Circular obstacle at (5,0), radius=1.2  (p_collision <= 0.005)
  - Box control chance constraints: |u_i| <= 0.15 per axis  (p_violation <= 0.005)
  - Covariance: P0 = diag([0.1,0.1,0.01,0.01]),  Pf = diag([0.04,0.04,0.01,0.01])

Run:
    python examples/example_constrained.py
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse

from sqrt_cs.systems.double_integrator import build_double_integrator
from sqrt_cs.problems.sqrt_qr_cs import SqrtQRCovarianceSteering
from sqrt_cs.utils.chance_constraints import AffineCC, CircularObstacle
from sqrt_cs.core.scp_params import SCPParams
from sqrt_cs.core.scvx_star import SCvxStar


def _make_init_ref(mu0, muf, N, dt, A, B, arc_height=2.5):
    """
    Build a dynamically-consistent initial mean/control reference that
    arcs above the obstacle (at x=5, y=arc_height).

    Trajectory: cosine arch so that vy=0 at both endpoints.
      py(k) = arc_height * (1 - cos(2π k/N)) / 2
      vy(k) = arc_height * π/T * sin(2π k/N)    (T = N*dt)
      uy(k) ≈ (vy(k+1) - vy(k)) / dt            (constant x-velocity → ux≈0)
    """
    T = N * dt
    k_arr = np.arange(N + 1)
    alpha = k_arr / N

    nx = len(mu0)
    nu = B.shape[1]

    init_mu = np.zeros((N + 1, nx))
    init_mu[:, 0] = alpha * muf[0]                                      # px linear
    init_mu[:, 1] = arc_height * (1 - np.cos(2 * np.pi * alpha)) / 2   # py cosine arch
    init_mu[:, 2] = muf[0] / T                                          # vx constant
    init_mu[:, 3] = arc_height * np.pi / T * np.sin(2 * np.pi * alpha) # vy

    # Control: v[k] = delta-velocity / dt  (approximately)
    init_v = np.zeros((N, nu))
    init_v[:, 0] = 0.0                                                   # ux ≈ 0
    init_v[:, 1] = np.diff(init_mu[:, 3]) / dt                          # uy

    return init_mu, init_v


def run(N: int = 20, verbose: bool = True):
    """
    Solve constrained CS for a 2-D double integrator with obstacle avoidance.

    Parameters
    ----------
    N       : horizon length  (MATLAB uses 20..160)
    verbose : print convergence table

    Returns
    -------
    problem : solved SqrtQRCovarianceSteering instance
    solver  : SCvxStar instance
    """
    # ---- System (matches MATLAB: nw=nx=4, G=sqrt(q*dt)*I) ----
    d = 2
    T_total = 30.0
    dt = T_total / N
    q = 0.005   # process noise spectral density (MATLAB single-obstacle scenario)

    A, B, _ = build_double_integrator(dt=dt, d=d)
    nx, nu = A.shape[0], B.shape[1]   # nx=4, nu=2

    nw = nx
    G = np.sqrt(q * dt) * np.eye(nw)   # (4,4), all-state noise

    # ---- Covariance boundary conditions (MATLAB) ----
    P0 = np.diag([0.1, 0.1, 0.01, 0.01])
    Pf = np.diag([0.04, 0.04, 0.01, 0.01])

    # ---- Mean boundary conditions ----
    mu0 = np.array([0.0, 0.0, 0.0, 0.0])
    muf = np.array([10.0, 0.0, 0.0, 0.0])

    # ---- Cost matrices ----
    Q = 0.1 * np.eye(nx)
    R = np.eye(nu)

    # ---- Chance constraints ----
    risk = 0.005
    u_max = 0.15

    # Control box: P(|u_x| <= u_max) >= 1-risk, P(|u_y| <= u_max) >= 1-risk
    # Affine CC: alpha' v[k] + z * ||L[k]' alpha||_2 <= beta
    cc_control = [
        AffineCC(alpha=np.array([ 1.0, 0.0]), beta=u_max, p=risk),
        AffineCC(alpha=np.array([-1.0, 0.0]), beta=u_max, p=risk),
        AffineCC(alpha=np.array([ 0.0, 1.0]), beta=u_max, p=risk),
        AffineCC(alpha=np.array([ 0.0,-1.0]), beta=u_max, p=risk),
    ]

    # Circular obstacle: center=(5,0), radius=1.2
    # pos_idx=[0,1] selects (px,py) from state [px,py,vx,vy]
    obstacle = CircularObstacle(
        center=np.array([5.0, 0.0]),
        radius=1.2,
        pos_idx=[0, 1],
        p=risk,
        nodes=list(range(1, N + 1)),
    )

    # ---- Dynamically-consistent initial reference ----
    init_mu, init_v = _make_init_ref(mu0, muf, N, dt, A, B, arc_height=2.5)

    # ---- Problem ----
    problem = SqrtQRCovarianceSteering(
        N=N,
        A=A, B=B, G=G,
        P0=P0, Pf=Pf,
        Q=Q, R=R,
        mu0=mu0, muf=muf,
        chance_constraints_control=cc_control,
        obstacles=[obstacle],
        init_mu=init_mu,
        init_v=init_v,
    )

    # ---- Solver params (matches MATLAB: tol_feas=1e-4, k_max=100) ----
    params = SCPParams(
        r_init=0.1,
        r_max=1.0,
        w_init=100.0,
        tol_opt=1e-3,
        tol_feas=1e-4,
        k_max=200,
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
        loss = problem.compute_covariance_propagation_loss(sol)
        print(f"Max cov propagation loss : {loss.max():.2e}")
        print(f"Mean cov propagation loss: {loss.mean():.2e}")

        # Check obstacle clearance (3-sigma ellipse from solved trajectory)
        if problem.mu_sol is not None:
            mu_arr = problem.mu_sol        # (N+1, nx)
            S_arr  = problem.S_sol         # (N+1, nx, nx)
            obs_center = obstacle.center
            obs_radius = obstacle.radius
            min_dist = np.inf
            for k in range(N + 1):
                pos = mu_arr[k, :2]
                dist = np.linalg.norm(pos - obs_center)
                min_dist = min(min_dist, dist)
            print(f"Min distance (mean) to obstacle center : {min_dist:.3f} m  (radius={obs_radius})")

    return problem, solver


def plot_results(problem: SqrtQRCovarianceSteering, solver: SCvxStar):
    """Plot mean trajectory, covariance ellipses, obstacle, and convergence."""
    it = solver.iter_state
    fig = plt.figure(figsize=(16, 5))
    ax_path = fig.add_subplot(131)
    ax_obj  = fig.add_subplot(132)
    ax_chi  = fig.add_subplot(133)

    # ---- Mean trajectory + covariance ellipses (pos subspace) ----
    if problem.mu_sol is not None and problem.P_sol is not None:
        mu  = problem.mu_sol    # (N+1, nx)
        P   = problem.P_sol     # (N+1, nx, nx)

        ax_path.plot(mu[:, 0], mu[:, 1], "b-o", markersize=3, label="Mean trajectory")
        ax_path.plot(mu[0, 0], mu[0, 1], "gs", markersize=8, label="Start")
        ax_path.plot(mu[-1, 0], mu[-1, 1], "r*", markersize=12, label="Goal")

        # 2-sigma covariance ellipses at each node (position subspace)
        z_2sigma = 2.0
        for k in range(problem.N + 1):
            Ppos = P[k][:2, :2]   # 2×2 position covariance
            vals, vecs = np.linalg.eigh(Ppos)
            angle = np.degrees(np.arctan2(vecs[1, -1], vecs[0, -1]))
            w, h = 2 * z_2sigma * np.sqrt(vals)
            ellipse = Ellipse(
                xy=(mu[k, 0], mu[k, 1]), width=w, height=h, angle=angle,
                edgecolor="cornflowerblue", facecolor="none", alpha=0.5, linewidth=0.8,
            )
            ax_path.add_patch(ellipse)

        # Obstacle
        obs_circle = plt.Circle((5, 0), 1.2, color="tomato", alpha=0.4, label="Obstacle")
        ax_path.add_patch(obs_circle)
        ax_path.add_patch(plt.Circle((5, 0), 1.2, color="tomato", fill=False, linewidth=2))

    ax_path.set_xlim(-0.5, 11)
    ax_path.set_ylim(-1.5, 4)
    ax_path.set_aspect("equal")
    ax_path.set_xlabel("x [m]")
    ax_path.set_ylabel("y [m]")
    ax_path.set_title("Mean trajectory + 2σ ellipses")
    ax_path.legend(fontsize=7)
    ax_path.grid(True, alpha=0.4)

    # ---- Objective convergence ----
    ax_obj.semilogy(it.hist_J0, color="steelblue")
    ax_obj.set_xlabel("SCP iteration")
    ax_obj.set_ylabel("Objective J0")
    ax_obj.set_title("Objective convergence")
    ax_obj.grid(True, alpha=0.4)

    # ---- Feasibility convergence ----
    ax_chi.semilogy(it.hist_chi, color="darkorange")
    ax_chi.axhline(1e-4, color="k", linestyle="--", linewidth=0.8, label="tol_feas")
    ax_chi.set_xlabel("SCP iteration")
    ax_chi.set_ylabel("Feasibility χ")
    ax_chi.set_title("Constraint violation")
    ax_chi.legend(fontsize=8)
    ax_chi.grid(True, alpha=0.4)

    plt.tight_layout()
    plt.savefig("constrained_results.png", dpi=150)
    print("Saved constrained_results.png")
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
