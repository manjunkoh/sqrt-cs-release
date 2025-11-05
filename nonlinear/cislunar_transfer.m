clc; clear;

addpath ../astrodynamics_base/
addpath(genpath("../cislunar-sequential-covariance-steering"))
addpath ../SCvxStar/src

addpath(genpath("./utils"))
addpath ./src

yalmip("clear")
figure_settings

load('../cislunar-sequential-covariance-steering/mat_files/NRHO_to_Halo_reference.mat', ...
    'x_ref0', 'u_ref0');

init_ref.x = x_ref0;
init_ref.u = u_ref0;

p = L2NRHOtoL1Halo_stochastic();

[A,B,~,G] = p.DS.discretize_LT_SDE(init_ref.x, init_ref.u, p.t_his, p.g);

% Extract dimensions
N = size(A, 3);  % Time horizon
nx = size(A, 1);
nu = size(B, 2);
nw = size(G, 2);

% Objective weights
Q = 0.01 * eye(nx);  % State cost matrix
R = 0.01 * eye(nu);   % Control cost matrix

%% Solve using FullCovarianceSteering
disp('=== Solving Cislunar Transfer with FullCovarianceSteering ===');

prob_fc = FullCovarianceSteering(...
    A=A, B=B, G=G, ...
    P_0=p.P0, P_f=p.P_fin, ...
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
disp('=== Solving Cislunar Transfer with SqrtQRCovarianceSteering ===');

% Scale the problem to make covariance matrices closer to identity
D = diag(1./sqrt(diag(p.P_fin)));
P0_s = D * p.P0 * D;
P_fin_s = D * p.P_fin * D;
A_sys_s = zeros(nx, nx, N);
B_sys_s = zeros(nx, nu, N);
G_sys_s = zeros(nx, nw, N);
for k = 1:N
    A_sys_s(:,:,k) = D * A(:,:,k) * inv(D);
    B_sys_s(:,:,k) = D * B(:,:,k);
    G_sys_s(:,:,k) = D * G(:,:,k);
end

% Initial guess for SqrtQRCovarianceSteering
% Interpolate Cholesky factors from initial to final covariance
init.S = interpolate_lower_triangular(chol(P0_s, 'lower'), chol(P_fin_s, 'lower'), N+1, 'log-cholesky');
init.L = zeros(nu, nx, N);

prob_qr = SqrtQRCovarianceSteering(init, ...
    N=N, ...
    A_sys=A_sys_s, B_sys=B_sys_s, G_sys=G_sys_s, ...
    P_0=P0_s, P_f=P_fin_s, Q=Q, R=R, ...
    objective_type='LQR');

% SCP parameters
scp_params = SCPParams();
scp_params.tol_opt = 1E-2;
scp_params.tol_feas = 1E-4;
scp_params.r_init = 0.1;
scp_params.penalty_method = 'ALwithL1';

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
        % solution_unscaled.K(:,:,k) = solution_unscaled.L(:,:,k) / solution_unscaled.S(:,:,k);
        solution_unscaled.K(:,:,k) = prob_qr.K(:,:,k) * D;
        % solution_unscaled.P_u(:,:,k) = solution_unscaled.L(:,:,k) * solution_unscaled.L(:,:,k)';
        solution_unscaled.P_u(:,:,k) = solution_unscaled.K(:,:,k) * solution_unscaled.P(:,:,k) * solution_unscaled.K(:,:,k)';
        solution_unscaled.L(:,:,k) = solution_unscaled.K(:,:,k) * solution_unscaled.S(:,:,k);
    end

    solution_unscaled.mu = prob_qr.sol.mu;
    solution_unscaled.v = prob_qr.sol.v;
    
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
ylabels = {'$x$, [LU]', '$y$, [LU]', '$z$, [LU]', '$v_x$, [VU]', '$v_y$, [VU]', '$v_z$, [VU]'};
figure
tiledlayout(nx, 1, "TileSpacing","tight", "Padding","tight");
for i = 1:nx
    nexttile
    plot(0:N, squeeze(sqrt(prob_fc.P(i,i,:))), 'k-', 'DisplayName', 'FullCovariance');
    hold on
    plot(0:N, squeeze(sqrt(solution_unscaled.P(i,i,:))), 'r-', 'DisplayName', 'SqrtQR');
    hold off

    if i == nx
        xlabel('Node $k$')
        legend('Location', 'best', 'Orientation','horizontal')
    end
    ylabel(ylabels{i})
end

exportgraphics(gcf, 'cislunar_transfer_covariance_diagonals.png')

%% Plot the control covariance diagonals for each method, for each control variable
figure
tiledlayout(nu, 1, "TileSpacing","tight", "Padding","tight");
for i = 1:nu
    nexttile
    plot(0:N-1, squeeze(sqrt(prob_fc.P_u(i,i,:))), 'k-', 'DisplayName', 'FullCovariance');
    hold on
    plot(0:N-1, squeeze(sqrt(solution_unscaled.P_u(i,i,:))), 'r-', 'DisplayName', 'SqrtQR');
    hold off
    if i == nu
        xlabel('Node $k$')
        legend('Location', 'best', 'Orientation','horizontal')
    end
    ylabel(ylabels{i})
    yscale log
end

% exportgraphics(gcf, 'cislunar_transfer_control_covariance_diagonals.png')