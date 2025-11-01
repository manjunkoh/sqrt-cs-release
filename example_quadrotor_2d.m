%% Example: Quadrotor 2D path planning
% Discrete triple integrator lateral/longitudinal dynamics
clear; clc;
addpath ../SCvxStar/src/
addpath(genpath('./utils'))
addpath('./src')

% Time horizon and discretization
N = 60;
dt = 0.1;

% Dimensions
nx = 6; % triple integrator in x and y (position and velocity and accel?)
nu = 2; % two inputs (x and y thrust)
nw = 6;

% Build block matrices as in the figure
I2 = eye(2);
Z2 = zeros(2);
A_block = [I2, dt*I2, 0*I2;
           Z2, I2,  dt*I2;
           Z2, Z2,  I2];
% A_sys is time-varying; replicate
A_sys = repmat(A_block, [1,1,N]);

B_block = [Z2; Z2; dt*I2];
B_sys = repmat(B_block, [1,1,N]);

G_block = 0.1 * eye(nx);
G_sys = repmat(G_block, [1,1,N]);

% Covariance boundary conditions
Sigma_i = eye(nx);
Sigma_f = 0.1 * eye(nx);

% Mean boundary conditions
mu_i = [20; 0; 0; 0; 0; 0];
mu_f = zeros(nx,1);

% Waypoint constraints: intermediate mean positions
% At node 20: position [13, 5], at node 40: position [7, -5]
% State is [x, y, vx, vy, ax, ay], so we only specify position (first 2 elements)
% Use NaN for unconstrained components (velocity and acceleration will be free)
waypoints = {
    struct('node', 20, 'mu', [13; 5; NaN; NaN; NaN; NaN]);  % Waypoint at node 20
    struct('node', 40, 'mu', [7; -5; NaN; NaN; NaN; NaN]);  % Waypoint at node 40
};

% Chance constraint parameters (boxes): example linear constraints
% State: [x, y, vx, vy, ax, ay] for 6D triple integrator
p_x = 0.005;
% State constraints: x in [-3, 25], y in [-7, 7]
alpha_x1 = [1; 0; 0; 0; 0; 0];  % x <= 25 (x is first element)
alpha_x2 = [-1; 0; 0; 0; 0; 0]; % x >= -3  ->  -x <= 3
alpha_x3 = [0; 1; 0; 0; 0; 0];  % y <= 7 (y is second element)
alpha_x4 = [0; -1; 0; 0; 0; 0]; % y >= -7  ->  -y <= 7
beta_x = [25, 3, 7, 7];

% Build chance constraint structs for BlockCholeskySteering
state_cc = {
    struct('type', 'affine', 'alpha', alpha_x1, 'beta', beta_x(1), 'p', p_x, 'nodes', 1:N+1);
    struct('type', 'affine', 'alpha', alpha_x2, 'beta', beta_x(2), 'p', p_x, 'nodes', 1:N+1);
    struct('type', 'affine', 'alpha', alpha_x3, 'beta', beta_x(3), 'p', p_x, 'nodes', 1:N+1);
    struct('type', 'affine', 'alpha', alpha_x4, 'beta', beta_x(4), 'p', p_x, 'nodes', 1:N+1);
};

% Objective weights
Q = 0.001 * eye(nx);
R = 0.01 * eye(nu);

%% Solve with BlockCholeskySteering
disp('=== Solving with BlockCholeskySteering ===');
prob_bc = BlockCholeskySteering(...
    A=A_sys, B=B_sys, G=G_sys, ...
    P_0=Sigma_i, P_f=Sigma_f, ...
    Q=Q, R=R, ...
    N=N, ...
    chance_constraints_state=state_cc, ...
    mu_0=mu_i, mu_f=mu_f, ...
    waypoints=waypoints);

tic
diagnostic_bc = prob_bc.solve('verbose', 0, 'solver', 'mosek');
time_bc = toc;

if diagnostic_bc.problem == 0 || diagnostic_bc.problem == 4
    P_bc = prob_bc.P;
    fprintf('BlockCholeskySteering solved successfully in %.3f seconds\n', time_bc);
    fprintf('  Objective value: %.6f\n', prob_bc.optimal_objective);
else
    disp('BlockCholeskySteering solver failed.');
    disp(yalmiperror(diagnostic_bc.problem));
end

%% Solve with FullCovarianceSteering
disp('=== Solving with FullCovarianceSteering ===');
% Reference covariance for linearization of chance constraints
P_ref = 1.2 * eye(nx);  % Time-invariant reference

prob_fc = FullCovarianceSteering(...
    A=A_sys, B=B_sys, G=G_sys, ...
    P_0=Sigma_i, P_f=Sigma_f, ...
    P_ref=P_ref, ...
    Q=Q, R=R, ...
    N=N, ...
    mu_0=mu_i, mu_f=mu_f, ...
    chance_constraints_state = state_cc, ...
    waypoints=waypoints);

tic
diagnostic_fc = prob_fc.solve('verbose', 0, 'solver', 'mosek');
time_fc = toc;

if diagnostic_fc.problem == 0 || diagnostic_fc.problem == 4
    P_fc = prob_fc.P;
    fprintf('FullCovarianceSteering solved successfully in %.3f seconds\n', time_fc);
    fprintf('  Objective value: %.6f\n', prob_fc.optimal_objective);
else
    disp('FullCovarianceSteering solver failed.');
    disp(yalmiperror(diagnostic_fc.problem));
    prob_fc = [];
end

%% Solve with SqrtQRCovarianceSteering
disp('=== Solving with SqrtQRCovarianceSteering ===');

% Initial guess for SqrtQRCovarianceSteering
init.S = interpolate_lower_triangular(chol(Sigma_i, 'lower'), chol(Sigma_f, 'lower'), N+1, 'log-cholesky');
init.L = zeros(nu, nx, N);
% Initialize mean trajectory to satisfy waypoints
init.mu = linspace_vec(mu_i, mu_f, N+1);
init.v = zeros(nu, N);

prob_qr = SqrtQRCovarianceSteering(init, ...
    N=N, ...
    A_sys=A_sys, B_sys=B_sys, G_sys=G_sys, ...
    P_0=Sigma_i, P_f=Sigma_f, Q=Q, R=R, ...
    objective_type='LQG', ...
    chance_constraints_state=state_cc, ...
    waypoints=waypoints, ...
    mu_0=mu_i, mu_f=mu_f);

scp_params = SCPParams();
scp_params.tol_opt = 1E-2;
scp_params.tol_feas = 1E-4;

tic
flag_solved_qr = prob_qr.solve(save_bool=false, scp_params=scp_params);
time_qr = seconds(prob_qr.scp.report.time);
J_qr = NaN;

if flag_solved_qr
    prob_qr.postprocess();
    % Compute objective value (LQG objective)
    J_qr = 0;
    for k = 1:prob_qr.N
        J_qr = J_qr + trace(Q * prob_qr.P(:,:,k)) + trace(R * prob_qr.P_u(:,:,k)) ...
            + prob_qr.mu(:,k)' * Q * prob_qr.mu(:,k) + prob_qr.v(:,k)' * R * prob_qr.v(:,k);
    end
    fprintf('SqrtQRCovarianceSteering solved successfully in %.3f seconds\n', time_qr);
    fprintf('  Number of iterations: %d\n', prob_qr.scp.report.iters);
    fprintf('  Objective value: %.6f\n', J_qr);
else
    disp('SqrtQRCovarianceSteering solver failed.');
    prob_qr = [];
end

%% Plot results
if (diagnostic_fc.problem == 0 || diagnostic_fc.problem == 4) || flag_solved_qr
    
    figure
    hold on
    
    % Plot constraint boundaries (box constraints: x in [-3, 25], y in [-7, 7])
    plot([-3, -3], [-10, 10], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Constraint: x = -3');
    plot([25, 25], [-10, 10], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Constraint: x = 25');
    plot([-5, 27], [-7, -7], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Constraint: y = -7');
    plot([-5, 27], [7, 7], 'r--', 'LineWidth', 1.5, 'DisplayName', 'Constraint: y = 7');
    
    % Extract position covariance (2x2 submatrix from 6x6 covariance)
    % State is [x, y, vx, vy, ax, ay], so position indices are 1:2
    idx_pos = 1:2;
    
    % Plot BlockCholeskySteering results
    % Plot covariance ellipses (every 5th step to avoid clutter)
    % for k = 1:5:N+1
    %     P_pos_bc = P_bc(idx_pos, idx_pos, k);
    %     plot3sigmaEllipse(prob_bc.mu(idx_pos, k), P_pos_bc, 'b', 'HandleVisibility', 'off');
    % end
    % % Plot mean trajectory
    % plot(prob_bc.mu(1,:), prob_bc.mu(2,:), 'b+-', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Block Cholesky');
    
    % Plot FullCovarianceSteering results
    if ~isempty(prob_fc) && (diagnostic_fc.problem == 0 || diagnostic_fc.problem == 4)
        % Plot covariance ellipses (every 5th step to avoid clutter)
        for k = 1:5:N+1
            P_pos_fc = P_fc(idx_pos, idx_pos, k);
            plot3sigmaEllipse(prob_fc.mu(idx_pos, k), P_pos_fc, 'g', 'HandleVisibility', 'off');
        end
        % Plot mean trajectory
        plot(prob_fc.mu(1,:), prob_fc.mu(2,:), 'g+-', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Full Covariance');
    end
    
    % Plot SqrtQRCovarianceSteering results
    if ~isempty(prob_qr) && flag_solved_qr
        % Plot covariance ellipses (every 5th step to avoid clutter)
        for k = 1:5:N+1
            plot3sigmaEllipse(prob_qr.mu(idx_pos, k), prob_qr.P(idx_pos, idx_pos, k), 'm', 'HandleVisibility', 'off');
        end
        % Plot mean trajectory
        plot(prob_qr.mu(1,:), prob_qr.mu(2,:), 'm+-', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Square Root QR');
    end
    
    % Plot initial and terminal conditions
    plot3sigmaEllipse(mu_i(idx_pos), Sigma_i(idx_pos, idx_pos), 'r', 'LineWidth', 2, 'DisplayName', 'Initial');
    plot3sigmaEllipse(mu_f(idx_pos), Sigma_f(idx_pos, idx_pos), 'r--', 'LineWidth', 2, 'DisplayName', 'Terminal');
    
    % Plot waypoints
    for i = 1:length(waypoints)
        wp = waypoints{i};
        plot(wp.mu(1), wp.mu(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', 'r', ...
            'DisplayName', sprintf('Waypoint (node %d)', wp.node));
    end
    
    xlabel('$x$ (position)', 'Interpreter', 'latex')
    ylabel('$y$ (position)', 'Interpreter', 'latex')
    title('Quadrotor 2D Path Planning')
    % legend('Location', 'best')
    grid on
    axis equal

end

%% Summary
fprintf('\n=== Solution Summary ===\n');
% if diagnostic_bc.problem == 0 || diagnostic_bc.problem == 4
%     fprintf('BlockCholeskySteering:  Time = %.3f s, Objective = %.6f\n', time_bc, prob_bc.optimal_objective);
% else
%     fprintf('BlockCholeskySteering:  Failed\n');
% end

if (diagnostic_fc.problem == 0 || diagnostic_fc.problem == 4) && ~isempty(prob_fc)
    fprintf('FullCovarianceSteering: Time = %.3f s, Objective = %.6f\n', time_fc, prob_fc.optimal_objective);
else
    fprintf('FullCovarianceSteering: Failed\n');
end

if flag_solved_qr && ~isempty(prob_qr)
    fprintf('SqrtQRCovarianceSteering: Time = %.3f s, Objective = %.6f, Iterations = %d\n', time_qr, J_qr, prob_qr.scp.report.iters);
else
    fprintf('SqrtQRCovarianceSteering: Failed\n');
end
