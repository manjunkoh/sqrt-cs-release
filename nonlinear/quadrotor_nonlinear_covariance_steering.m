% NonlinearSqrtQRCovarianceSteering with Quadrotor2D system
% Based on parameters from the paper

clc; clear;

addpath ../SCvxStar/src/
addpath ../astrodynamics_base/
addpath(fullfile(fileparts(mfilename('fullpath')), 'src'))
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')))

yalmip("clear")
figure_settings

%% System Parameters (from image/paper)
% Quadrotor2D system: state is [x, y, vx, vy] (4D)
nx = 4;  % State dimension: position (x,y) and velocity (vx, vy)
nu = 2;  % Control dimension: thrust in x and y directions
nw = 2;  % Process noise dimension

% Time horizon and discretization
N = 25;           % Number of discrete steps (from image)
sigma = 15;       % Time scale (from image)
t_final = sigma;  % Final time
dt = t_final / N; % Time step
t_his = linspace(0, t_final, N+1);

% Drag coefficient
c_d = 0.005;  % k_D from image

% Create Quadrotor2D dynamical system
DS = Quadrotor2D(c_d);

%% Covariance Steering Problem Parameters (from image)
% Initial state distribution
mu_0 = [1; 8; 2; 0];        % Initial mean: [x, y, vx, vy]
P0 = 0.01 * eye(nx);        % Initial covariance: 0.01 * I

% Terminal state distribution
mu_f = [1; 2; -1; 0];       % Terminal mean: [x, y, vx, vy]
P_fin = 0.1 * eye(nx);      % Terminal covariance: 0.1 * I

% Process noise
gamma = 0.01;  % Noise scale from image
g = @(t,x) gamma * [zeros(2); eye(2)];  % Process noise input matrix

% Cost matrices (from image)
% Mean cost: l(x_tau, u_tau) = 10 * ||u_tau||^2
% Q_x,tau = 5 * I, Q_u,tau = I
Q = 5 * eye(nx);    % State cost matrix
R = 10 * eye(nu);   % Control cost matrix (10 * I for mean cost 10 * ||u||^2)

% Terminal mean error weight (from image)
w_xf = 1000;  % Not directly used in standard formulation, but noted

%% Chance Constraint (from image Equation 71)
% P(||e1 * xi_tau||_1 <= 6) >= 0.9
% where e1 = [1, 0] (selects first component)
% This becomes: P(|x| <= 6) >= 0.9, which is P(-6 <= x <= 6) >= 0.9
% We can express this as two affine chance constraints:
% P(x <= 6) >= 0.9 and P(-x <= 6) >= 0.9 (i.e., P(x >= -6) >= 0.9)
e1 = [1; 0; 0; 0];  % Selects x position (first state component)
p_chance = 0.1;     % Violation probability: 1 - 0.9 = 0.1
beta_constraint = 6;  % Constraint value

% State chance constraint: |x| <= 6 with probability >= 0.9
state_chance_constraint = {
    struct('type', 'affine', 'alpha', e1, 'beta', beta_constraint, 'p', p_chance, 'nodes', 1:N+1),  % x <= 6
    struct('type', 'affine', 'alpha', -e1, 'beta', beta_constraint, 'p', p_chance, 'nodes', 1:N+1), % x >= -6 (i.e., -x <= 6)
};

%% Initial Guess for NonlinearSqrtQRCovarianceSteering
% Initial control guess from image: u_hat_k^1 = [-0.3 -0.1]^T
v_guess = repmat([-0.3; -0.1], 1, N);  % Constant initial control guess

% Initial state trajectory guess (linear interpolation between initial and final)
x_guess = linspace_vec(mu_0, mu_f, N+1);

% Initial guess for covariance square roots
init_guess = struct(...
    'S', interpolate_lower_triangular(chol(P0, 'lower'), chol(P_fin, 'lower'), N+1, 'log-cholesky'), ...
    'L', zeros(nu, nx, N), ...
    'mu', x_guess, ...
    'v', v_guess);

%% Solve using NonlinearSqrtQRCovarianceSteering
disp('=== Solving Quadrotor2D Transfer with NonlinearSqrtQRCovarianceSteering ===');

prob_nl = NonlinearSqrtQRCovarianceSteering(init_guess, ...
    N=N, ...
    DS=DS, ...
    g=g, ...
    t_his=t_his, ...
    P_0=P0, P_f=P_fin, ...
    objective_type='LQG', Q=Q, R=R, ...
    mu_0=mu_0, mu_f=mu_f, ...
    chance_constraints_state=state_chance_constraint);

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

% Plot 1: Mean position trajectories
subplot(2,3,1);
hold on;
plot(0:N, x_guess(1,:), 'k--', 'DisplayName', 'Initial Guess (x)', 'LineWidth', 1);
plot(0:N, x_guess(2,:), 'k:', 'DisplayName', 'Initial Guess (y)', 'LineWidth', 1);
plot(0:N, prob_nl.mu(1,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Nonlinear x Mean');
plot(0:N, prob_nl.mu(2,:), 'g-', 'LineWidth', 2, 'DisplayName', 'Nonlinear y Mean');
xlabel('Time Step $k$');
ylabel('Position');
title('Mean Position Trajectories');
legend('Location', 'best');
grid on;

% Plot 2: Mean velocity trajectories
subplot(2,3,2);
hold on;
plot(0:N, x_guess(3,:), 'k--', 'DisplayName', 'Initial Guess (vx)', 'LineWidth', 1);
plot(0:N, x_guess(4,:), 'k:', 'DisplayName', 'Initial Guess (vy)', 'LineWidth', 1);
plot(0:N, prob_nl.mu(3,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Nonlinear vx Mean');
plot(0:N, prob_nl.mu(4,:), 'g-', 'LineWidth', 2, 'DisplayName', 'Nonlinear vy Mean');
xlabel('Time Step $k$');
ylabel('Velocity');
title('Mean Velocity Trajectories');
legend('Location', 'best');
grid on;

% Plot 3: Mean control trajectory
subplot(2,3,3);
hold on;
plot(0:N-1, v_guess(1,:), 'k--', 'LineWidth', 1, 'DisplayName', 'Initial Guess');
plot(0:N-1, prob_nl.v(1,:), 'b-', 'LineWidth', 2, 'DisplayName', 'u_x');
plot(0:N-1, prob_nl.v(2,:), 'r-', 'LineWidth', 2, 'DisplayName', 'u_y');
xlabel('Time Step $k$');
ylabel('Control Value');
title('Mean Control Trajectory');
legend('Location', 'best');
grid on;

% Plot 4: State covariance diagonals (position)
subplot(2,3,4);
hold on;
plot(0:N, squeeze(sqrt(prob_nl.P(1,1,:))), 'b-', ...
    'LineWidth', 1.5, 'DisplayName', '$\sigma_x$');
plot(0:N, squeeze(sqrt(prob_nl.P(2,2,:))), 'g-', ...
    'LineWidth', 1.5, 'DisplayName', '$\sigma_y$');
xlabel('Time Step $k$');
ylabel('Standard Deviation');
title('Position Covariance (Diagonal)');
legend('Location', 'best');
grid on;

% Plot 5: State covariance diagonals (velocity)
subplot(2,3,5);
hold on;
plot(0:N, squeeze(sqrt(prob_nl.P(3,3,:))), 'b-', ...
    'LineWidth', 1.5, 'DisplayName', '$\sigma_{vx}$');
plot(0:N, squeeze(sqrt(prob_nl.P(4,4,:))), 'g-', ...
    'LineWidth', 1.5, 'DisplayName', '$\sigma_{vy}$');
xlabel('Time Step $k$');
ylabel('Standard Deviation');
title('Velocity Covariance (Diagonal)');
legend('Location', 'best');
grid on;

% Plot 6: Control covariance
subplot(2,3,6);
hold on;
plot(0:N-1, squeeze(sqrt(prob_nl.P_u(1,1,:))), 'b-', ...
    'LineWidth', 1.5, 'DisplayName', '$\sigma_{u_x}$');
plot(0:N-1, squeeze(sqrt(prob_nl.P_u(2,2,:))), 'r-', ...
    'LineWidth', 1.5, 'DisplayName', '$\sigma_{u_y}$');
xlabel('Time Step $k$');
ylabel('Standard Deviation');
title('Control Covariance');
legend('Location', 'best');
grid on;
yscale log;

sgtitle(sprintf('Quadrotor2D NonlinearSqrtQRCovarianceSteering (Converged: %d iterations)', ...
    prob_nl.scp.report.iters));

% exportgraphics(gcf, 'figures/quadrotor_nonlinear_covariance_steering.png')

%% Plot 2D trajectory with uncertainty ellipses
figure;
hold on;
% Plot trajectory
plot(prob_nl.mu(1,:), prob_nl.mu(2,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Mean Trajectory');

% Plot uncertainty ellipses at selected nodes
ellipse_nodes = 1:5:N+1;  % Every 5th node
for k = ellipse_nodes
    mu_k = prob_nl.mu([1,2], k);  % Position only
    P_k = prob_nl.P([1,2], [1,2], k);  % Position covariance
    plotEllipse(mu_k, P_k, 3, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);  % 3-sigma ellipse
end

% Plot initial and final positions
plot(prob_nl.mu(1,1), prob_nl.mu(2,1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot(prob_nl.mu(1,end), prob_nl.mu(2,end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'End');

% Plot chance constraint boundary (x = ±6)
xline(6, 'r--', 'LineWidth', 1, 'DisplayName', 'Chance Constraint ($x = 6$)');
xline(-6, 'r--', 'LineWidth', 1, 'DisplayName', 'Chance Constraint ($x = -6$)');

xlabel('$x$ position');
ylabel('$y$ position');
title('2D Trajectory with Uncertainty Ellipses (3$\sigma$)');
legend('Location', 'best');
grid on;
axis equal;

exportgraphics(gcf, 'figures/quadrotor_nonlinear_covariance_steering_trajectory.png')

fprintf('\n=== Script Complete ===\n');

