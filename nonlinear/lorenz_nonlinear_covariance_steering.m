% NonlinearSqrtQRCovarianceSteering with Lorenz system

clc; clear;

addpath ../SCvxStar/src/
addpath ../astrodynamics_base/
addpath(fullfile(fileparts(mfilename('fullpath')), 'example4'))
addpath(fullfile(fileparts(mfilename('fullpath')), 'src'))
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')))

yalmip("clear")
figure_settings

%% System Parameters
% Load initial reference trajectory from ExampleClass4
disp('=== Computing Initial Reference Trajectory ===');
prob_ref = ExampleClass4Optimizer();
SCvxParams = SCPParams();
SCvxParams.tol_opt = 1E-3;
prob_ref.solve(scp_params=SCvxParams);

% Extract system dimensions
nx = prob_ref.DS.nx;  % 3 for Lorenz system
nu = prob_ref.DS.nu;  % 1 for Lorenz system
N = prob_ref.Nseg;    % Number of segments
t_his = prob_ref.t_his;
DeltaT = prob_ref.DeltaT;

%% Define process noise
nw = nx;  % Process noise dimension
process_noise_scale = 0.01;
g = @(t,x) process_noise_scale * eye(nw);  % Process noise input matrix

%% Covariance Steering Problem Parameters
P0 = 0.1 * eye(nx);      % Initial covariance
P_fin = 0.05 * eye(nx);  % Final covariance
Q = 0.1 * eye(nx);       % State cost matrix
R = 0.1 * eye(nu);       % Control cost matrix

% Initial and final mean states
mu_0 = prob_ref.x_init;
mu_f = prob_ref.x_fin;

%% Initial Guess for NonlinearSqrtQRCovarianceSteering
% Use reference trajectory from ExampleClass4 as initial guess
x_guess = prob_ref.sol.x;  % Initial mean state trajectory guess
v_guess = prob_ref.sol.u;  % Initial mean control trajectory guess

% Initial guess for covariance square roots
init_guess = struct(...
    'S', interpolate_lower_triangular(chol(P0, 'lower'), chol(P_fin, 'lower'), N+1, 'log-cholesky'), ...
    'L', zeros(nu, nx, N), ...
    'mu', x_guess, ...
    'v', v_guess);

%% Chance Constraints
control_chance_constraint = {struct('type', 'norm', 'gamma', ExampleClass4.u_max, 'p', 0.1, 'n', nu)};

%% Solve using NonlinearSqrtQRCovarianceSteering
disp('=== Solving Lorenz Transfer with NonlinearSqrtQRCovarianceSteering ===');

prob_nl = NonlinearSqrtQRCovarianceSteering(init_guess, ...
    N=N, ...
    DS=prob_ref.DS, ...
    g=g, ...
    t_his=t_his, ...
    P_0=P0, P_f=P_fin, ...
    objective_type='DV99', Q=Q, R=R,...
    mu_0=mu_0, mu_f=mu_f,...
    chance_constraints_control=control_chance_constraint);

% SCP parameters
scp_params = SCPParams();
scp_params.tol_opt = 1E-3;
scp_params.tol_feas = 1E-4;
scp_params.linearization = 'inexact';
% scp_params.penalty_method = 'ALwithL1';

tic
flag_solved = prob_nl.solve(save_bool=false, scp_params=scp_params);
time_total = seconds(prob_nl.scp.report.time);
J_final = NaN;

%% Post-processing and Results
if flag_solved
    prob_nl.postprocess();
    
    % Compute final objective value
    J_final = 0;
    for k = 1:N
        J_final = J_final + trace(Q * prob_nl.P(:,:,k)) + trace(R * prob_nl.P_u(:,:,k));
        J_final = J_final + prob_nl.mu(:,k)' * Q * prob_nl.mu(:,k) + prob_nl.v(:,k)' * R * prob_nl.v(:,k);
    end
    J_final = J_final + prob_nl.mu(:,N+1)' * Q * prob_nl.mu(:,N+1);
    
    fprintf('NonlinearSqrtQRCovarianceSteering solved successfully in %.3f seconds\n', time_total);
    fprintf('  Objective value: %.6f\n', J_final);
    fprintf('  Number of SCP iterations: %d\n', prob_nl.scp.report.iters);
    fprintf('Solution complete. Results stored in prob_nl object.\n');
else
    disp('NonlinearSqrtQRCovarianceSteering solver failed.');
    prob_nl = [];
    return
end

%% Plotting
figure();

% Plot 1: Mean state trajectories
subplot(2,2,1);
hold on;
plot(0:N, prob_ref.sol.x(1,:), 'k--', 'DisplayName', 'Initial Ref (X)', 'LineWidth', 1);
plot(0:N, prob_ref.sol.x(2,:), 'k:', 'DisplayName', 'Initial Ref (Y)', 'LineWidth', 1);
plot(0:N, prob_ref.sol.x(3,:), 'k-.', 'DisplayName', 'Initial Ref (Z)', 'LineWidth', 1);
plot(0:N, prob_nl.mu(1,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Nonlinear X Mean');
plot(0:N, prob_nl.mu(2,:), 'g-', 'LineWidth', 2, 'DisplayName', 'Nonlinear Y Mean');
plot(0:N, prob_nl.mu(3,:), 'r-', 'LineWidth', 2, 'DisplayName', 'Nonlinear Z Mean');
xlabel('Time Step $k$');
ylabel('State Value');
title('Mean State Trajectories');
legend('Location', 'best');
grid on;

% Plot 2: Mean control trajectory
subplot(2,2,2);
hold on;
plot(0:N-1, prob_ref.sol.u(1,:), 'k--', 'LineWidth', 1, 'DisplayName', 'Initial Ref');
plot(0:N-1, prob_nl.v(1,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Nonlinear Mean Control');
xlabel('Time Step $k$');
ylabel('Control Value');
title('Mean Control Trajectory');
legend('Location', 'best');
grid on;

% Plot 3: State covariance diagonals
subplot(2,2,3);
hold on;
ylabels = {'$x$', '$y$', '$z$'};
colors = {'b', 'g', 'r'};
for j = 1:nx
    plot(0:N, squeeze(sqrt(prob_nl.P(j,j,:))), [colors{j} '-'], ...
        'LineWidth', 1.5, 'DisplayName', sprintf('$\\sigma_{%s}$', ylabels{j}));
end
xlabel('Time Step $k$');
ylabel('Standard Deviation');
title('State Covariance (Diagonal)');
legend('Location', 'best');
grid on;

% Plot 4: Control covariance
subplot(2,2,4);
hold on;
plot(0:N-1, squeeze(sqrt(prob_nl.P_u(1,1,:))), 'b-', ...
    'LineWidth', 1.5, 'DisplayName', 'Control $\\sigma$');
xlabel('Time Step $k$');
ylabel('Standard Deviation');
title('Control Covariance');
legend('Location', 'best');
grid on;
yscale log;

sgtitle(sprintf('NonlinearSqrtQRCovarianceSteering (Converged: %d iterations)', ...
    prob_nl.scp.report.iters));

% exportgraphics(gcf, 'figures/lorenz_nonlinear_covariance_steering.png')

%% Plot 3D trajectory
figure;
view(3);
hold on;
plot3d(prob_ref.sol.x, 'k--', 'DisplayName', 'Initial Reference', 'LineWidth', 1.5);
plot3d(prob_nl.mu', 'b-', 'DisplayName', 'Nonlinear Mean Trajectory', 'LineWidth', 2);
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
title('3D Trajectory Comparison');
legend('Location', 'best');
axis equal;
grid on;
view(45, 20);

exportgraphics(gcf, 'figures/lorenz_nonlinear_covariance_steering_3d.png')

fprintf('\n=== Script Complete ===\n');

