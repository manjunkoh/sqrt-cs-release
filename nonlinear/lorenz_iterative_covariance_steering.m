% Iterative Covariance Steering with Lorenz system

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
process_noise_scale = 0.1;
g = @(t,x) process_noise_scale * eye(nw);  % Process noise input matrix

%% Covariance Steering Problem Parameters
P0 = 0.1 * eye(nx);      % Initial covariance
P_fin = 0.05 * eye(nx);  % Final covariance
Q = 0.1 * eye(nx);       % State cost matrix
R = 0.1 * eye(nu);       % Control cost matrix

%% Iterative Covariance Steering (iCS) Parameters
i_max = 100;         % Maximum number of outer iterations
tol = 1e-3;         % Convergence tolerance for mean control
solver_method = 'SqrtQRCovarianceSteering';
% solver_method = 'FullCovarianceSteering';

%% Initial Guess for Outer Loop (û¹_k, K¹_k)
% Use reference trajectory from ExampleClass4 as initial mean control guess
v_hat_i = prob_ref.sol.u;  % û¹_k (mean control trajectory)
x_hat_i = prob_ref.sol.x;  % Initial mean state trajectory guess
K_hat_i = zeros(nu, nx, N);  % K¹_k (initial feedback gains - zero for now)

% Initial state is fixed
mu_0 = prob_ref.x_init;
mu_f = prob_ref.x_fin;

%% Outer Loop: Iterative Covariance Steering (iCS)
disp('=== Starting Iterative Covariance Steering (iCS) ===');
converged = false;
time_iCS_start = tic;

% Initialize warm-start structure for inner covariance steering problem
% After first iteration, this will be updated with the previous solution
inner_init = struct('S', interpolate_lower_triangular(chol(P0, 'lower'), chol(P_fin, 'lower'), N+1, 'log-cholesky'), ...
    'L', zeros(nu, nx, N), ...
    'mu', x_hat_i, ...
    'v', v_hat_i);

% SCP parameters for inner loop, used for SqrtQR
scp_params_inner = SCPParams();
scp_params_inner.tol_opt = 1E-3;
scp_params_inner.tol_feas = 1E-4;
scp_params_inner.r_init = 0.1;

control_chance_constraint = {struct('type', 'norm', 'gamma', ExampleClass4.u_max, 'p', 0.01, 'n', nu)};

for i = 1:i_max
    fprintf('\n--- Outer Iteration %d/%d ---\n', i, i_max);
    
    % Step 2: Propagate nonlinear mean dynamics with ûⁱ_k, Kⁱ_k
    % This generates x̄_k from the nonlinear system
    fprintf('  Propagating nonlinear mean dynamics...\n');
    % Propagate using nonlinear dynamics with mean control v_hat_i
    u_handle = @(t) getZOH(t, t_his, v_hat_i);
    [~, x_bar_k] = prob_ref.DS.propagate_with_LT(mu_0, t_his, u_handle);
    x_bar_k = x_bar_k';
    
    % Step 3: x̂ⁱ_k ← x̄_k
    x_hat_i = x_bar_k;
    fprintf('  Mean state trajectory updated.\n');
    
    % Step 4: Linearize about (x̂ⁱ, ûⁱ) & Step 5: Discretize
    fprintf('  Linearizing and discretizing system...\n');
    [A, B, ~, G] = prob_ref.DS.discretize_LT_SDE(x_hat_i, v_hat_i, t_his, g);
    
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
                objective_type='DV99', ...
                mu_0=mu_0, mu_f=mu_f, ...
                chance_constraints_control=control_chance_constraint, ...
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
                mu_0=mu_0, mu_f=mu_f);
            
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
        % Extract numerical values from sdpvar objects
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
            % If mean variables weren't defined, we need to propagate mean dynamics
            % using the new feedback gains. For now, use previous guess as approximation.
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
subplot(2,2,1);
hold on;
plot(0:N, prob_ref.sol.x(1,:), 'k--', 'DisplayName', 'Initial Ref (X)', 'LineWidth', 1);
plot(0:N, prob_ref.sol.x(2,:), 'k:', 'DisplayName', 'Initial Ref (Y)', 'LineWidth', 1);
plot(0:N, prob_ref.sol.x(3,:), 'k-.', 'DisplayName', 'Initial Ref (Z)', 'LineWidth', 1);

if strcmp(solver_method, 'SqrtQRCovarianceSteering') && isfield(prob_cs, 'mu') && ~isempty(prob_cs.mu)
    plot(0:N, prob_cs.mu(1,:), 'b-', 'LineWidth', 2, 'DisplayName', 'iCS X Mean');
    plot(0:N, prob_cs.mu(2,:), 'g-', 'LineWidth', 2, 'DisplayName', 'iCS Y Mean');
    plot(0:N, prob_cs.mu(3,:), 'r-', 'LineWidth', 2, 'DisplayName', 'iCS Z Mean');
else
    plot(0:N, x_hat_i(1,:), 'b-', 'LineWidth', 2, 'DisplayName', 'iCS X Mean');
    plot(0:N, x_hat_i(2,:), 'g-', 'LineWidth', 2, 'DisplayName', 'iCS Y Mean');
    plot(0:N, x_hat_i(3,:), 'r-', 'LineWidth', 2, 'DisplayName', 'iCS Z Mean');
end

xlabel('Time Step $k$');
ylabel('State Value');
title('Mean State Trajectories');
legend('Location', 'best');
grid on;

% Plot 2: Mean control trajectory
subplot(2,2,2);
hold on;
plot(0:N-1, prob_ref.sol.u(1,:), 'k--', 'LineWidth', 1, 'DisplayName', 'Initial Ref');
plot(0:N-1, v_hat_i(1,:), 'b-', 'LineWidth', 2, 'DisplayName', 'iCS Mean Control');
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
    plot(0:N, squeeze(sqrt(prob_cs.P(j,j,:))), [colors{j} '-'], ...
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
plot(0:N-1, squeeze(sqrt(prob_cs.P_u(1,1,:))), 'b-', ...
    'LineWidth', 1.5, 'DisplayName', 'Control $\\sigma$');
xlabel('Time Step $k$');
ylabel('Standard Deviation');
title('Control Covariance');
legend('Location', 'best');
grid on;
yscale log;

sgtitle(sprintf('iCS - %s (Converged: %d iterations)', ...
    solver_method, i));

% exportgraphics(gcf, 'figures/lorenz_iCS.png')

%% Plot 3D trajectory
figure;
view(3);
hold on;
plot3d(prob_ref.sol.x, 'k--', 'DisplayName', 'Initial Reference', 'LineWidth', 1.5);
if strcmp(solver_method, 'SqrtQRCovarianceSteering') && isfield(prob_cs, 'mu') && ~isempty(prob_cs.mu)
    plot3d(prob_cs.mu', 'b-', 'DisplayName', 'iCS Mean Trajectory', 'LineWidth', 2);
else
    plot3d(x_hat_i', 'b-', 'DisplayName', 'iCS Mean Trajectory', 'LineWidth', 2);
end
xlabel('$x$');
ylabel('$y$');
zlabel('$z$');
title('3D Trajectory Comparison');
legend('Location', 'best');
axis equal;
grid on;
view(45, 20);

exportgraphics(gcf, 'figures/lorenz_iCS_3d.png')

fprintf('\n=== Script Complete ===\n');

