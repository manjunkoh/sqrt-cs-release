% Iterative Covariance Steering with Quadrotor2D system
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
nw = 2; % Process noise dimension

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

%% Iterative Covariance Steering (iCS) Parameters
i_max = 100;         % Maximum number of outer iterations
tol = 1e-3;         % Convergence tolerance for mean control
solver_method = 'SqrtQRCovarianceSteering';
% solver_method = 'FullCovarianceSteering';

%% Initial Guess for Outer Loop
% Initial control guess from image: u_hat_k^1 = [-0.3 -0.1]^T
v_hat_i = repmat([-0.3; -0.1], 1, N);  % Constant initial control guess

% Initial state trajectory guess (linear interpolation between initial and final)
x_hat_i = linspace_vec(mu_0, mu_f, N+1);
K_hat_i = zeros(nu, nx, N);  % Initial feedback gains - zero for now

%% Outer Loop: Iterative Covariance Steering (iCS)
disp('=== Starting Iterative Covariance Steering (iCS) for Quadrotor2D ===');
converged = false;
time_iCS_start = tic;

% Initialize warm-start structure for inner covariance steering problem
inner_init = struct('S', interpolate_lower_triangular(chol(P0, 'lower'), chol(P_fin, 'lower'), N+1, 'log-cholesky'), ...
    'L', zeros(nu, nx, N), ...
    'mu', x_hat_i, ...
    'v', v_hat_i);

% SCP parameters for inner loop
scp_params_inner = SCPParams();
scp_params_inner.tol_opt = 1E-3;
scp_params_inner.tol_feas = 1E-4;
scp_params_inner.r_init = 0.1;

for i = 1:i_max
    fprintf('\n--- Outer Iteration %d/%d ---\n', i, i_max);
    
    % Step 2: Propagate nonlinear mean dynamics with ûⁱ_k, Kⁱ_k
    % This generates x̄_k from the nonlinear system
    fprintf('  Propagating nonlinear mean dynamics...\n');
    % Propagate using nonlinear dynamics with mean control v_hat_i
    u_handle = @(t) getZOH(t, t_his, v_hat_i);
    [~, x_bar_k] = DS.propagate_with_LT(mu_0, t_his, u_handle);
    x_bar_k = x_bar_k';
    
    % Step 3: x̂ⁱ_k ← x̄_k
    x_hat_i = x_bar_k;
    fprintf('  Mean state trajectory updated.\n');
    
    % Step 4: Linearize about (x̂ⁱ, ûⁱ) & Step 5: Discretize
    fprintf('  Linearizing and discretizing system...\n');
    [A, B, ~, G] = DS.discretize_LT_SDE(x_hat_i, v_hat_i, t_his, g);
    
    % Step 6: Solve Covariance Steering Problem
    fprintf('  Solving covariance steering problem (inner loop)...\n');
    
    % Solve using chosen method
    switch solver_method

        case 'SqrtQRCovarianceSteering'

            % Warm-start: Update mean trajectory guess with current reference
            inner_init.mu = x_hat_i;
            inner_init.v = v_hat_i;
            % S and L are warm-started from previous iteration (or initial guess for first iteration)

            if i == 1
                impose_mean_trust_region = false;
            else
                impose_mean_trust_region = true;
            end
            prob_cs = SqrtQRCovarianceSteering(inner_init, ...
                N=N, ...
                A_sys=A, B_sys=B, G_sys=G, ...
                P_0=P0, P_f=P_fin, Q=Q, R=R, ...
                objective_type='LQG', ...
                mu_0=mu_0, mu_f=mu_f, ...
                chance_constraints_state=state_chance_constraint, ...
                impose_mean_trust_region=impose_mean_trust_region, ...
                mean_trust_region_radius=0.1, ...
                mu_ref = x_bar_k);
            
            flag_solved = prob_cs.solve(scp_params=scp_params_inner);
            inner_solve_time = seconds(prob_cs.scp.report.time);
        
        case 'FullCovarianceSteering'
            prob_cs = FullCovarianceSteering(...
                A=A, B=B, G=G, ...
                P_0=P0, P_f=P_fin, ...
                Q=Q, R=R, ...
                N=N, ...
                mu_0=mu_0, mu_f=mu_f, ...
                chance_constraints_state=state_chance_constraint);
            
            inner_solve_start = tic;
            diagnostic = prob_cs.solve('verbose', 0, 'solver', 'mosek');
            inner_solve_time = toc(inner_solve_start);
            flag_solved = (diagnostic.problem == 0 || diagnostic.problem == 4);
            
        otherwise
            error('Invalid solver_method: %s', solver_method);
    end
    
    if ~flag_solved
        warning('Inner covariance steering problem failed at iteration %d.', i);
        break;
    end
    
    % Post-process to extract results
    if strcmp(solver_method, 'SqrtQRCovarianceSteering')
        prob_cs.postprocess();
        v_star = prob_cs.v;  % ū*_k (optimal mean control)
        K_star = prob_cs.K;  % K* (optimal feedback gains)
        
        % Warm-start update: Store solution for next iteration
        inner_init.S = value(prob_cs.sol.S);  % Warm-start S from current solution
        inner_init.L = value(prob_cs.sol.L);  % Warm-start L from current solution
        inner_init.mu = prob_cs.mu;           % Warm-start mean trajectory
        inner_init.v = prob_cs.v;             % Warm-start mean control
    else
        % For FullCovarianceSteering, extract mean control and feedback gains
        prob_cs.set_feedback_gains();
        K_star = prob_cs.K;
        
        % Extract mean control v if it exists, otherwise use previous guess
        if isfield(prob_cs, 'v') && ~isempty(prob_cs.v)
            v_star = prob_cs.v;
        else
            warning('FullCovarianceSteering: Mean control not available. Using previous guess.');
            v_star = v_hat_i;
        end
    end
    
    fprintf('  Inner problem solved in %.3f s.\n', inner_solve_time);
    
    % Step 8: Convergence Check
    % Check: max_k ||ū*_k - ûⁱ_k|| ≤ tol
    delta_v = max(vecnorm(v_star - v_hat_i, 2, 1));
    fprintf('  Delta V (max control change): %.6e (tol: %.6e)\n', delta_v, tol);
    
    if delta_v <= tol
        fprintf('  ✓ Convergence achieved!\n');
        converged = true;
        break;
    end
    
    % Step 11: Update reference for next iteration
    % ûⁱ⁺¹_k ← ū*_k, Kⁱ⁺¹_k ← K*
    v_hat_i = v_star;
    K_hat_i = K_star;
end

time_iCS_total = toc(time_iCS_start);

%% Post-processing and Results
fprintf('\n=== iCS Solution Summary ===\n');
fprintf('Total iCS solve time: %.3f seconds\n', time_iCS_total);
fprintf('Number of outer iterations: %d\n', i);
if converged
    fprintf('✓ Converged successfully!\n');
else
    fprintf('⚠ Did not converge.\n');
    return
end

% Compute final objective value
if strcmp(solver_method, 'SqrtQRCovarianceSteering')
    J_final = 0;
    for k = 1:N
        J_final = J_final + trace(Q * prob_cs.P(:,:,k)) + trace(R * prob_cs.P_u(:,:,k));
        if isfield(prob_cs, 'mu') && ~isempty(prob_cs.mu)
            J_final = J_final + prob_cs.mu(:,k)' * Q * prob_cs.mu(:,k) + prob_cs.v(:,k)' * R * prob_cs.v(:,k);
        end
    end
    if isfield(prob_cs, 'mu') && ~isempty(prob_cs.mu)
        J_final = J_final + prob_cs.mu(:,N+1)' * Q * prob_cs.mu(:,N+1);
    end
else
    J_final = prob_cs.optimal_objective;
end
fprintf('Final Objective Value: %.6f\n', J_final);

%% Plotting
figure();

% Plot 1: Mean state trajectories
subplot(2,3,1);
hold on;
plot(0:N, x_hat_i(1,:), 'k--', 'DisplayName', 'Initial Ref (x)', 'LineWidth', 1);
plot(0:N, x_hat_i(2,:), 'k:', 'DisplayName', 'Initial Ref (y)', 'LineWidth', 1);

if strcmp(solver_method, 'SqrtQRCovarianceSteering') && isfield(prob_cs, 'mu') && ~isempty(prob_cs.mu)
    plot(0:N, prob_cs.mu(1,:), 'b-', 'LineWidth', 2, 'DisplayName', 'iCS x Mean');
    plot(0:N, prob_cs.mu(2,:), 'g-', 'LineWidth', 2, 'DisplayName', 'iCS y Mean');
else
    plot(0:N, x_hat_i(1,:), 'b-', 'LineWidth', 2, 'DisplayName', 'iCS x Mean');
    plot(0:N, x_hat_i(2,:), 'g-', 'LineWidth', 2, 'DisplayName', 'iCS y Mean');
end
xlabel('Time Step $k$');
ylabel('Position');
title('Mean Position Trajectories');
legend('Location', 'best');
grid on;

% Plot 2: Mean velocity trajectories
subplot(2,3,2);
hold on;
plot(0:N, x_hat_i(3,:), 'k--', 'DisplayName', 'Initial Ref (vx)', 'LineWidth', 1);
plot(0:N, x_hat_i(4,:), 'k:', 'DisplayName', 'Initial Ref (vy)', 'LineWidth', 1);

if strcmp(solver_method, 'SqrtQRCovarianceSteering') && isfield(prob_cs, 'mu') && ~isempty(prob_cs.mu)
    plot(0:N, prob_cs.mu(3,:), 'b-', 'LineWidth', 2, 'DisplayName', 'iCS vx Mean');
    plot(0:N, prob_cs.mu(4,:), 'g-', 'LineWidth', 2, 'DisplayName', 'iCS vy Mean');
else
    plot(0:N, x_hat_i(3,:), 'b-', 'LineWidth', 2, 'DisplayName', 'iCS vx Mean');
    plot(0:N, x_hat_i(4,:), 'g-', 'LineWidth', 2, 'DisplayName', 'iCS vy Mean');
end
xlabel('Time Step $k$');
ylabel('Velocity');
title('Mean Velocity Trajectories');
legend('Location', 'best');
grid on;

% Plot 3: Mean control trajectory
subplot(2,3,3);
hold on;
plot(0:N-1, v_hat_i(1,:), 'b-', 'LineWidth', 2, 'DisplayName', 'u_x');
plot(0:N-1, v_hat_i(2,:), 'r-', 'LineWidth', 2, 'DisplayName', 'u_y');
xlabel('Time Step $k$');
ylabel('Control Value');
title('Mean Control Trajectory');
legend('Location', 'best');
grid on;

% Plot 4: State covariance diagonals (position)
subplot(2,3,4);
hold on;
if strcmp(solver_method, 'SqrtQRCovarianceSteering')
    plot(0:N, squeeze(sqrt(prob_cs.P(1,1,:))), 'b-', ...
        'LineWidth', 1.5, 'DisplayName', '$\sigma_x$');
    plot(0:N, squeeze(sqrt(prob_cs.P(2,2,:))), 'g-', ...
        'LineWidth', 1.5, 'DisplayName', '$\sigma_y$');
end
xlabel('Time Step $k$');
ylabel('Standard Deviation');
title('Position Covariance (Diagonal)');
legend('Location', 'best');
grid on;

% Plot 5: State covariance diagonals (velocity)
subplot(2,3,5);
hold on;
if strcmp(solver_method, 'SqrtQRCovarianceSteering')
    plot(0:N, squeeze(sqrt(prob_cs.P(3,3,:))), 'b-', ...
        'LineWidth', 1.5, 'DisplayName', '$\sigma_{vx}$');
    plot(0:N, squeeze(sqrt(prob_cs.P(4,4,:))), 'g-', ...
        'LineWidth', 1.5, 'DisplayName', '$\sigma_{vy}$');
end
xlabel('Time Step $k$');
ylabel('Standard Deviation');
title('Velocity Covariance (Diagonal)');
legend('Location', 'best');
grid on;

% Plot 6: 2D trajectory
subplot(2,3,6);
hold on;
plot(x_hat_i(1,:), x_hat_i(2,:), 'k--', 'DisplayName', 'Initial Reference', 'LineWidth', 1.5);
if strcmp(solver_method, 'SqrtQRCovarianceSteering') && isfield(prob_cs, 'mu') && ~isempty(prob_cs.mu)
    plot(prob_cs.mu(1,:), prob_cs.mu(2,:), 'b-', 'DisplayName', 'iCS Mean Trajectory', 'LineWidth', 2);
    % Plot initial and final positions
    plot(prob_cs.mu(1,1), prob_cs.mu(2,1), 'go', 'MarkerSize', 10, 'DisplayName', 'Start');
    plot(prob_cs.mu(1,end), prob_cs.mu(2,end), 'ro', 'MarkerSize', 10, 'DisplayName', 'End');
else
    plot(x_hat_i(1,:), x_hat_i(2,:), 'b-', 'DisplayName', 'iCS Mean Trajectory', 'LineWidth', 2);
end
xlabel('$x$ position');
ylabel('$y$ position');
title('2D Trajectory');
legend('Location', 'best');
grid on;
axis equal;

sgtitle(sprintf('Quadrotor2D iCS - %s (Converged: %d iterations)', ...
    solver_method, i));

% exportgraphics(gcf, 'figures/quadrotor_iCS.png')

%% Plot 2D trajectory with uncertainty ellipses
figure;
hold on;
if strcmp(solver_method, 'SqrtQRCovarianceSteering') && isfield(prob_cs, 'mu') && ~isempty(prob_cs.mu)
    % Plot trajectory
    plot(prob_cs.mu(1,:), prob_cs.mu(2,:), 'b-', 'LineWidth', 2, 'DisplayName', 'Mean Trajectory');
    
    % Plot uncertainty ellipses at selected nodes
    ellipse_nodes = 1:5:N+1;  % Every 5th node
    for k = ellipse_nodes
        mu_k = prob_cs.mu([1,2], k);  % Position only
        P_k = prob_cs.P([1,2], [1,2], k);  % Position covariance
        plotEllipse(mu_k, P_k, 3, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);  % 3-sigma ellipse
    end
    
    % Plot initial and final positions
    plot(prob_cs.mu(1,1), prob_cs.mu(2,1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
    plot(prob_cs.mu(1,end), prob_cs.mu(2,end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'End');
    
    % Plot chance constraint boundary (x = ±6)
    xline(6, 'r--', 'LineWidth', 1, 'DisplayName', 'Chance Constraint ($x = 6$)');
    xline(-6, 'r--', 'LineWidth', 1, 'DisplayName', 'Chance Constraint ($x = -6$)');
end
xlabel('$x$ position');
ylabel('$y$ position');
title('2D Trajectory with Uncertainty Ellipses (3$\sigma$)');
legend('Location', 'best');
grid on;
axis equal;

exportgraphics(gcf, 'figures/quadrotor_iCS_trajectory.png')

fprintf('\n=== Script Complete ===\n');

