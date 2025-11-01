%% Example: Cone-shaped path planning problem
% Solves with both BlockCholeskySteering (Okamoto et al. 2019) and 
% SqrtQRCovarianceSteering methods
clear; 
clc;
addpath ../SCvxStar/src/
addpath(genpath('./utils'))
addpath('./src')
yalmip('clear')

figure_settings

% Problem parameters
dt = 0.2;
N = 20;

% State: x = [x, y, vx, vy]'
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
Q = diag([1e-2, 1e-2, 1e-3, 1e-3]);
R = diag([1, 1]);

% Cone-shaped chance constraints: -0.2*(x-1) <= y <= 0.2*(x-1)
% Rearranged to alpha' * x <= beta form for state = [x, y, vx, vy]':
% Constraint 1: y <= 0.2(x-1)   ->  y <= 0.2x - 0.2  ->  0.2x + y <= 0.2
% Constraint 2: y >= -0.2(x-1)  ->  y >= -0.2x + 0.2  ->  0.2x - y <= 0.2
alpha1 = [0.2; 1; 0; 0]; beta1 = 0.2;  % Upper bound: y <= 0.2(x-1)
alpha2 = [0.2; -1; 0; 0]; beta2 = 0.2; % Lower bound: y >= -0.2(x-1)

% Violation probability for chance constraints
p_violation = 0.005;

% Build chance constraint structs (apply to all time steps k=1:N+1)
state_cc = {
    struct('type', 'affine', 'alpha', alpha1, 'beta', beta1, 'p', p_violation, 'nodes', 1:N+1);
    struct('type', 'affine', 'alpha', alpha2, 'beta', beta2, 'p', p_violation, 'nodes', 1:N+1);
};

%% Solve with BlockCholeskySteering (Okamoto 2019)
disp('=== Solving with BlockCholeskySteering (Okamoto 2019) ===');
prob_bc = BlockCholeskySteering(...
    A=A_sys, B=B_sys, G=G_sys, ...
    P_0=Sigma0, P_f=SigmaN, ...
    Q=Q, R=R, ...
    N=N, ...
    chance_constraints_state=state_cc, ...
    mu_0=mu0, mu_f=muN);

tic
diagnostic_bc = prob_bc.solve('verbose', 0, 'solver', 'mosek');
time_bc = toc;

if diagnostic_bc.problem == 0 || diagnostic_bc.problem == 4
    P_bc = prob_bc.P;
    fprintf('BlockCholeskySteering solved successfully in %.3f seconds\n', time_bc);
else
    disp('BlockCholeskySteering solver failed.');
    disp(yalmiperror(diagnostic_bc.problem));
    prob_bc = [];
end

%% Solve with SqrtQRCovarianceSteering
disp('=== Solving with SqrtQRCovarianceSteering ===');
% Initial guess for SqrtQRCovarianceSteering
init.S = interpolate_lower_triangular(chol(Sigma0, 'lower'), chol(SigmaN, 'lower'), N+1, 'log-cholesky');
init.L = zeros(nu, nx, N);
init.mu = zeros(nx, N+1);
init.v = zeros(nu, N);

prob_qr = SqrtQRCovarianceSteering(init, ...
    N=N, ...
    A_sys=A_sys, B_sys=B_sys, G_sys=G_sys, ...
    P_0=Sigma0, P_f=SigmaN, Q=Q, R=R, ...
    objective_type='LQR', ...
    chance_constraints_state=state_cc, ...
    mu_0=mu0, mu_f=muN);

scp_params = SCPParams();
scp_params.tol_opt = 1E-4;
scp_params.tol_feas = 1E-4;
scp_params.r_init = 1.0;

tic
flag_solved_qr = prob_qr.solve(save_bool=false, scp_params=scp_params);
time_qr = toc;

if flag_solved_qr
    prob_qr.postprocess();
    fprintf('SqrtQRCovarianceSteering solved successfully in %.3f seconds\n', time_qr);
    fprintf('  Number of iterations: %d\n', prob_qr.scp.report.iters);
else
    disp('SqrtQRCovarianceSteering solver failed.');
    prob_qr = [];
end

%% Plot results
figure
hold on

% Plot constraint boundaries (cone walls)
% State is [x, vx, y, vy]', so constraint alpha'*state <= beta becomes:
% alpha(1)*x + alpha(2)*vx + alpha(3)*y + alpha(4)*vy <= beta
% For plotting in (x,y) space, we set vx=0, vy=0 and solve for y:
% alpha(1)*x + alpha(2)*y = beta -> y = (beta - alpha(1)*x) / alpha(2) if alpha(2) != 0
x_ch = linspace(-12, 1, 100)';
for j = 1:length(state_cc)
    cc = state_cc{j};
    if strcmp(cc.type, 'affine')
        a = cc.alpha;
        b = cc.beta;
        y_ch = (b - a(1)*x_ch) / a(2);
        plot(x_ch, y_ch, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Constraint boundary');
    end
end

% Plot results from BlockCholeskySteering
if ~isempty(prob_bc)
    % Plot covariance ellipses
    for k = 1:N+1
        plot3sigmaEllipse(prob_bc.mu(:,k), P_bc(:,:,k), 'b', 'HandleVisibility', 'off')
    end
    % Plot mean trajectory
    plot(prob_bc.mu(1,:), prob_bc.mu(2,:), 'b+-', 'LineWidth', 1.5, 'DisplayName', 'Block Cholesky');
end

% Plot results from SqrtQRCovarianceSteering
if ~isempty(prob_qr)
    % Plot covariance ellipses
    for k = 1:N+1
        plot3sigmaEllipse(prob_qr.mu(:,k), prob_qr.P(:,:,k), 'g', 'HandleVisibility', 'off')
    end
    % Plot mean trajectory
    plot(prob_qr.mu(1,:), prob_qr.mu(2,:), 'g+-', 'LineWidth', 1.5, 'DisplayName', 'Square Root QR');
end

% Plot initial and terminal constraints
plot3sigmaEllipse(mu0, Sigma0, 'r', 'LineWidth', 2, 'DisplayName', 'Initial');
plot3sigmaEllipse(muN, SigmaN, 'r--', 'LineWidth', 2, 'DisplayName', 'Terminal');

axis equal
xlabel("$x_1$", 'Interpreter', 'latex')
ylabel("$x_2$", 'Interpreter', 'latex')
title('Covariance Steering with Cone-Shaped Chance Constraints')
legend('Location', 'best')
grid on
% exportgraphics(gcf, './figures/example_cone_path.png')
