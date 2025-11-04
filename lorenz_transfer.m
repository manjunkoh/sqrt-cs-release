clc; clear;

addpath(fullfile(fileparts(mfilename('fullpath')), 'example4'))
addpath(fullfile(fileparts(mfilename('fullpath')), 'src'))
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'utils')))

yalmip("clear")
figure_settings

%% Load or compute reference trajectory from ExampleClass4
disp('=== Computing Reference Trajectory ===');
prob_ref = ExampleClass4Optimizer();
SCvxParams = SCPParams();
SCvxParams.tol_opt = 1E-3;
prob_ref.solve(scp_params=SCvxParams);

% Extract reference trajectory
init_ref.x = prob_ref.sol.x;
init_ref.u = prob_ref.sol.u;

% System dimensions
nx = prob_ref.DS.nx;  % 3 for Lorenz system
nu = prob_ref.DS.nu;  % 1 for Lorenz system
N = prob_ref.Nseg;    % Number of segments

%% Define process noise
% Process noise input matrix: add noise to all states
nw = nx;  % Process noise dimension
g = eye(nw);  % Process noise input matrix (identity - noise affects all states)
% Scale process noise by a factor (adjust as needed)
process_noise_scale = 0.1;
g = @(t,x) process_noise_scale * g;

%% Discretize system with process noise
[A, B, ~, G] = prob_ref.DS.discretize_LT_SDE(init_ref.x, init_ref.u, prob_ref.t_his, g);

%% Define covariance boundary conditions
% Initial covariance (small, well-defined)
P0 = 0.1 * eye(nx);

% Final covariance (smaller than initial, indicating reduction in uncertainty)
P_fin = 0.05 * eye(nx);

%% Objective weights
Q = 0.1 * eye(nx);  % State cost matrix
R = 0.1 * eye(nu);   % Control cost matrix

%% Solve using FullCovarianceSteering
disp('=== Solving Lorenz Transfer with FullCovarianceSteering ===');

prob_fc = FullCovarianceSteering(...
    A=A, B=B, G=G, ...
    P_0=P0, P_f=P_fin, ...
    Q=Q, R=R, ...
    N=N);

tic
diagnostic_fc = prob_fc.solve('verbose', 1, 'solver', 'mosek');
time_fc = toc;

if diagnostic_fc.problem == 0 || diagnostic_fc.problem == 4
    fprintf('FullCovarianceSteering solved successfully in %.3f seconds\n', time_fc);
    fprintf('  Objective value: %.6f\n', prob_fc.optimal_objective);
    fprintf('Solution complete. Results stored in prob_fc object.\n');
else
    disp('FullCovarianceSteering solver failed.');
    disp(yalmiperror(diagnostic_fc.problem));
    prob_fc = [];
end

%% Solve using SqrtQRCovarianceSteering
disp('=== Solving Lorenz Transfer with SqrtQRCovarianceSteering ===');

% Scale the problem to make covariance matrices closer to identity
D = diag(1./sqrt(diag(P_fin)));
P0_s = D * P0 * D;
P_fin_s = D * P_fin * D;
A_sys_s = zeros(nx, nx, N);
B_sys_s = zeros(nx, nu, N);
G_sys_s = zeros(nx, nw, N);
for k = 1:N
    A_sys_s(:,:,k) = D * A(:,:,k) * inv(D);
    B_sys_s(:,:,k) = D * B(:,:,k);
    G_sys_s(:,:,k) = D * G(:,:,k);
end
Q_s = inv(D) * Q * inv(D);
R_s = R;

% Initial guess for SqrtQRCovarianceSteering
% Interpolate Cholesky factors from initial to final covariance
init.S = interpolate_lower_triangular(chol(P0_s, 'lower'), chol(P_fin_s, 'lower'), N+1, 'log-cholesky');
init.L = zeros(nu, nx, N);

% Use reference trajectory as initial guess for mean (if needed)
% For pure covariance steering, we don't need mean variables
init.mu = zeros(nx, N+1);
init.v = zeros(nu, N);

prob_qr = SqrtQRCovarianceSteering(init, ...
    N=N, ...
    A_sys=A_sys_s, B_sys=B_sys_s, G_sys=G_sys_s, ...
    P_0=P0_s, P_f=P_fin_s, Q=Q_s, R=R_s, ...
    objective_type='LQG');

% SCP parameters
scp_params = SCPParams();
scp_params.tol_opt = 1E-2;
scp_params.tol_feas = 1E-4;
scp_params.r_init = 0.1;
% scp_params.penalty_method = 'ALwithL1';

tic
flag_solved_qr = prob_qr.solve(save_bool=false, scp_params=scp_params);
time_qr = seconds(prob_qr.scp.report.time);
J_qr = NaN;

%%
if flag_solved_qr
    prob_qr.postprocess();

    % Unscale the solution
    solution_unscaled = struct();
    for k = 1:prob_qr.N+1
        solution_unscaled.S(:,:,k) = inv(D) * prob_qr.sol.S(:,:,k);
        solution_unscaled.P(:,:,k) = solution_unscaled.S(:,:,k) * solution_unscaled.S(:,:,k)';
    end
    for k = 1:prob_qr.N
        solution_unscaled.K(:,:,k) = prob_qr.K(:,:,k) * D;
        solution_unscaled.P_u(:,:,k) = solution_unscaled.K(:,:,k) * solution_unscaled.P(:,:,k) * solution_unscaled.K(:,:,k)';
        solution_unscaled.L(:,:,k) = solution_unscaled.K(:,:,k) * solution_unscaled.S(:,:,k);
    end

    % Mean variables (if they exist)
    if isfield(prob_qr.sol, 'mu') && ~isempty(prob_qr.sol.mu)
        solution_unscaled.mu = prob_qr.sol.mu;
        solution_unscaled.v = prob_qr.sol.v;
    else
        solution_unscaled.mu = zeros(nx, N+1);
        solution_unscaled.v = zeros(nu, N);
    end
    
    % Compute objective value
    J_qr = prob_qr.objective(solution_unscaled);
    
    fprintf('SqrtQRCovarianceSteering solved successfully in %.3f seconds\n', time_qr);
    fprintf('  Objective value: %.6f\n', J_qr);
    fprintf('  Number of SCP iterations: %d\n', prob_qr.scp.report.iters);
    fprintf('Solution complete. Results stored in prob_qr object.\n');
else
    disp('SqrtQRCovarianceSteering solver failed.');
    prob_qr = [];
end

%% Summary
fprintf('\n=== Solution Summary ===\n');
if ~isempty(prob_fc) && (diagnostic_fc.problem == 0 || diagnostic_fc.problem == 4)
    fprintf('FullCovarianceSteering:  Time = %.3f s, Objective = %.6f\n', time_fc, prob_fc.optimal_objective);
else
    fprintf('FullCovarianceSteering:  Failed\n');
end

if ~isempty(prob_qr) && flag_solved_qr
    fprintf('SqrtQRCovarianceSteering:  Time = %.3f s, Objective = %.6f, Iterations = %d\n', ...
        time_qr, J_qr, prob_qr.scp.report.iters);
else
    fprintf('SqrtQRCovarianceSteering:  Failed\n');
end

%% Plot the covariance diagonals for each method, for each state variable
ylabels = {'$x$', '$y$', '$z$'};
figure
tiledlayout(nx, 1, "TileSpacing","tight", "Padding","tight");
for i = 1:nx
    nexttile
    plot(0:N, squeeze(sqrt(prob_fc.P(i,i,:))), 'k-', 'DisplayName', 'FullCovariance', 'LineWidth', 1.5);
    hold on
    if ~isempty(prob_qr) && flag_solved_qr
        plot(0:N, squeeze(sqrt(solution_unscaled.P(i,i,:))), 'r--', 'DisplayName', 'SqrtQR', 'LineWidth', 1.5);
    end
    hold off
    grid on

    if i == nx
        xlabel('Node $k$')
        legend('Location', 'best', 'Orientation','horizontal')
    end
    ylabel(ylabels{i})
    title(sprintf('State Covariance: %s', ylabels{i}))
end

exportgraphics(gcf, 'figures/lorenz_transfer_covariance_diagonals.png')

%% Plot the control covariance diagonals for each method, for each control variable
figure
tiledlayout(nu, 1, "TileSpacing","tight", "Padding","tight");
for i = 1:nu
    nexttile
    plot(0:N-1, squeeze(sqrt(prob_fc.P_u(i,i,:))), 'k-', 'DisplayName', 'FullCovariance', 'LineWidth', 1.5);
    hold on
    if ~isempty(prob_qr) && flag_solved_qr
        plot(0:N-1, squeeze(sqrt(solution_unscaled.P_u(i,i,:))), 'r--', 'DisplayName', 'SqrtQR', 'LineWidth', 1.5);
    end
    hold off
    grid on
    if i == nu
        xlabel('Node $k$')
        legend('Location', 'best', 'Orientation','horizontal')
    end
    ylabel('Control Covariance')
    title('Control Covariance')
    yscale log
end

exportgraphics(gcf, 'figures/lorenz_transfer_control_covariance_diagonals.png')

%% Plot the reference trajectory in 3D
figure;
view(3)
hold on;
plot3d(init_ref.x, 'b--', DisplayName='Reference trajectory', LineWidth=1.5)
xlabel('$x$')
ylabel('$y$')
zlabel('$z$')
title('Lorenz System Reference Trajectory')
legend()
axis equal
grid on
view(45, 20)

exportgraphics(gcf, 'figures/lorenz_transfer_reference_trajectory.png')

%%
prob_qr.scp.plot_time()
prob_qr.scp.print_time_info()