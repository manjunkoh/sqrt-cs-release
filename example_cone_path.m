%% Example: Cone-shaped path planning problem (from figure)
% Implements the system and constraints described in the attached figure.
clear; clc;
addpath ../SCvxStar/src/
figure_settings

% Problem parameters
dt = 0.2;
N = 20;

% State: x = [x, vx, y, vy]'
nx = 4;
nu = 2;
nw = 4; % process noise dimension same as state here

% Build A and B (time-invariant)
A = [1  0   dt  0;
     0  1   0   dt;
     0  0   1   0;
     0  0   0   1];

B = [0.5*dt^2  0;
     0         0.5*dt^2;
     dt        0;
     0         dt];

% Replicate for time-varying arrays
A_sys = repmat(A, [1,1,N]);
B_sys = repmat(B, [1,1,N]);

% Process noise (diagonal)
G = diag([0.01, 0.01, 0.01, 0.01]);
G_sys = repmat(G, [1,1,N]);

% Covariance boundary conditions
Sigma0 = diag([0.1, 0.1, 0.01, 0.01]);
SigmaN = 0.5 * Sigma0;

% Mean boundary conditions
mu0 = [-10; 1; 0; 0];
muN = zeros(nx,1);

% Objective weights
Q = diag([10, 10, 1, 1]);
R = diag([1e3, 1e3]);

% Cone-shaped linear constraints: 0.2*(x-1) <= y <= -0.2*(x-1)
% Rearranged to alpha' * x <= beta form:
% y <= -0.2 x + 0.2  ->  [0.2, 0, 1, 0]' * state <= 0.2
% y >= 0.2 x - 0.2   ->  [0.2, 0, -1, 0]' * state <= 0.2
alpha1 = [0.2; -1; 0; 0]; beta1 = 0.2; % y <= -0.2(x-1)
alpha2 = [0.2; 1; 0; 0]; beta2 = 0.2; % y >=  0.2(x-1)

% Probabilities for chance constraints
p_state = 0.005;

% Build chance constraint lists
% state_cc = cell();
% apply cone constraints at all time steps
state_cc = [
    AffineChanceConstraint(alpha1, beta1, p_state, 1:N);
    AffineChanceConstraint(alpha2, beta2, p_state, 1:N);
];

% Initial guess
init.S = interpolate_lower_triangular(chol(Sigma0, 'lower'), chol(SigmaN, 'lower'), N+1, 'log-cholesky');
init.L = zeros(nu, nx, N);
init.mu = zeros(nx, N+1);
init.v = zeros(nu, N);

% Options
prob = SqrtQRCovarianceSteering(init, ...
    N=N, nx=nx, nu=nu, nw=nw, ...
    A_sys=A_sys, ...
    B_sys=B_sys, ...
    G_sys=G_sys, ...
    P_0=Sigma0, P_f=SigmaN, Q=Q, R=R, ...
    objective_type='LQG', ...
    chance_constraints_state=state_cc, ...
    mu_0=mu0, mu_f=muN);

flag_solved = prob.solve();

if flag_solved
    prob.postprocess();
    disp('Solved successfully.');
else
    disp('Solver failed.');
end
prob.postprocess();
%%
figure
hold on
x_ch = linspace(-10, 1)';

for j = 1:length(state_cc)
    a = state_cc(j).a;
    b = state_cc(j).b;
    y_ch = (b - a(1)*x_ch) / a(2);
    line(x_ch,y_ch,'Color','r', 'LineWidth', 1)
end

for k = 1:N
    plot3sigmaEllipse(prob.mu(:,k), prob.P(:,:,k), 'b')
end
plot3sigmaEllipse(mu0, Sigma0, 'r')
plot3sigmaEllipse(muN, SigmaN, 'r--')
axis equal
xlabel("$x_1$")
ylabel("$x_2$")
exportgraphics(gcf, 'example_cone_path.png')

%%
prob.scp.plot_iter_history()
exportgraphics(gcf, 'example_cone_path_iter_history.png')
%%
prob.scp.plot_time()
exportgraphics(gcf, 'example_cone_path_time.png')