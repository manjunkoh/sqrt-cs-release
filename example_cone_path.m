%% Constrained covariance steering from Okamoto 2018
% Solves with both BlockCholeskySteering (Okamoto et al. 2019) and 
% SqrtQRCovarianceSteering methods
clc; clear;

addpath ./SCvxStar/src/
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

sdp_settings = sdpsettings('verbose', 0, 'solver', 'mosek');

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
diagnostic_bc = prob_bc.solve(sdp_settings);
time_bc = toc;

if diagnostic_bc.problem == 0 || diagnostic_bc.problem == 4
    P_bc = prob_bc.P;
    fprintf('BlockCholeskySteering solved successfully in %.3f seconds\n', time_bc);
else
    disp('BlockCholeskySteering solver failed.');
    disp(yalmiperror(diagnostic_bc.problem));
    prob_bc = [];
end

%% Solve with FullCovarianceSteering (iterative wrapper)
disp('=== Solving with FullCovarianceSteering (iterative) ===');
% The wrapper automatically:
% 1. Solves without chance constraints to get initial P_ref, Y_ref
% 2. Iteratively solves with chance constraints, updating references

wrapper_fc = FullCovarianceSteeringIterative(...
    A=A_sys, B=B_sys, G=G_sys, ...
    P_0=Sigma0, P_f=SigmaN, ...
    Q=Q, R=R, ...
    N=N, ...
    chance_constraints_state=state_cc, ...
    mu_0=mu0, mu_f=muN, ...
    max_iters=20, ...
    tol_opt=1e-4, ...
    tol_feas=1e-4, ...
    verbose=true);

tic
[diagnostic_fc, prob_fc] = wrapper_fc.solve(sdp_settings);
time_fc = toc;

if diagnostic_fc.problem == 0 || diagnostic_fc.problem == 4
    P_fc = prob_fc.P;
    fprintf('FullCovarianceSteering solved successfully in %.3f seconds\n', time_fc);
    fprintf('  Number of iterations: %d\n', length(wrapper_fc.iter_history) - 1);
else
    disp('FullCovarianceSteering solver failed.');
    disp(yalmiperror(diagnostic_fc.problem));
    prob_fc = [];
end

%% Solve with SqrtQRCovarianceSteering
disp('=== Solving with SqrtQRCovarianceSteering ===');

scp_params = SCPParams();
scp_params.tol_opt = 1E-4;
scp_params.tol_feas = 1E-4;

% Initial guess for SqrtQRCovarianceSteering
clear init

init.S = interpolate_lower_triangular(chol(Sigma0, 'lower'), chol(SigmaN, 'lower'), N+1, 'log-cholesky');
init.L = zeros(nu, nx, N);
init.mu = linspace_vec(mu0, muN, N+1);
init.v = zeros(nu, N);

init.t_L = zeros(1, N);
init.t_S = zeros(1, N);
init.t_mu = zeros(1, N);
init.t_v = zeros(1, N);

for k = 1:N
	init.t_L(k) = trace(init.L(:,:,k) * init.L(:,:,k)' * R);
	init.t_S(k) = trace(init.S(:,:,k) * init.S(:,:,k)' * Q);
	init.t_mu(k) = init.mu(:,k)' * Q * init.mu(:,k);
	init.t_v(k) = init.v(:,k)' * R * init.v(:,k);
end

prob_qr = SqrtQRCovarianceSteeringOptimizer(init, ...
    N=N, ...
    A_sys=A_sys, B_sys=B_sys, G_sys=G_sys, ...
    P_0=Sigma0, P_f=SigmaN, Q=Q, R=R, ...
    objective_type='LQR', ...
    chance_constraints_state=state_cc, ...
    mu_0=mu0, mu_f=muN);

flag_solved_qr = prob_qr.solve(save_bool=false, scp_params=scp_params);
time_qr = seconds(prob_qr.scp.report.time);

if flag_solved_qr
    prob_qr.postprocess();
    fprintf('SqrtQRCovarianceSteering solved successfully in %.3f seconds\n', time_qr);
    fprintf('  Number of iterations: %d\n', prob_qr.scp.report.iters);
else
    disp('SqrtQRCovarianceSteering solver failed.');
    prob_qr = [];
end

%% Plot results
figure(Position=[0, 0, 15, 15])
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

% Plot results from FullCovarianceSteering
if ~isempty(prob_fc)
    % Plot covariance ellipses
    for k = 1:N+1
        plot3sigmaEllipse(prob_fc.mu(:,k), P_fc(:,:,k), 'm', 'HandleVisibility', 'off')
    end
    % Plot mean trajectory
    plot(prob_fc.mu(1,:), prob_fc.mu(2,:), 'm+-', 'LineWidth', 1.5, 'DisplayName', 'Full Covariance');
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
% plot3sigmaEllipse(mu0, Sigma0, 'r', 'LineWidth', 2, 'DisplayName', 'Initial');
plot3sigmaEllipse(muN, SigmaN, 'r--', 'LineWidth', 2, 'DisplayName', 'Terminal');

axis equal
xlabel("$x_1$", 'Interpreter', 'latex')
ylabel("$x_2$", 'Interpreter', 'latex')
legend(legendUnq(), 'Location', 'best')
grid on
% exportgraphics(gcf, './figures/example_cone_path.png')

%% Comparison of results
fprintf('\n=== Comparison of Results ===\n');
fprintf('%-40s %-25s %-25s %-25s\n', 'Metric', 'BlockCholeskySteering', 'FullCovarianceSteering', 'SqrtQRCovarianceSteering');
fprintf('%s\n', repmat('-', 1, 115));

% Initialize values
time_bc_str = 'Failed';
time_fc_str = 'Failed';
time_qr_str = 'Failed';
iters_fc_str = 'N/A';
iters_qr_str = 'N/A';
J_bc_str = 'Failed';
J_fc_str = 'Failed';
J_qr_str = 'Failed';
err_term_bc_str = 'Failed';
err_term_fc_str = 'Failed';
err_term_qr_str = 'Failed';
err_mu_bc_str = 'Failed';
err_mu_fc_str = 'Failed';
err_mu_qr_str = 'Failed';
max_viol_bc_str = 'Failed';
max_viol_fc_str = 'Failed';
max_viol_qr_str = 'Failed';

% Compute BlockCholeskySteering metrics
if ~isempty(prob_bc)
    time_bc_str = sprintf('%.4f', time_bc);
    J_bc = prob_bc.optimal_objective;
    J_bc_str = sprintf('%.6f', J_bc);
    
    P_term_bc = prob_bc.P(:,:,end);
    % Compute error as largest eigenvalue of (P_term_bc - SigmaN)
    % If largest eigenvalue <= 0, achieved covariance is smaller/equal to target, error = 0
    diff_cov_bc = P_term_bc - SigmaN;
    diff_cov_bc = (diff_cov_bc + diff_cov_bc') / 2; % Ensure symmetry
    max_eig_bc = max(eig(diff_cov_bc));
    err_term_bc = max(0, max_eig_bc);  % Error is 0 if max_eig <= 0
    err_term_bc_str = sprintf('%.6e', err_term_bc);
    
    mu_term_bc = prob_bc.mu(:,end);
    err_mu_bc = norm(mu_term_bc - muN);
    err_mu_bc_str = sprintf('%.6e', err_mu_bc);
    
    % Chance constraint satisfaction
    inv_norm_p = norminv(1 - p_violation);
    max_viol_bc = 0;
    for j = 1:length(state_cc)
        cc = state_cc{j};
        if strcmp(cc.type, 'affine')
            alpha = cc.alpha;
            beta = cc.beta;
            for k = 1:N+1
                mu_k = prob_bc.mu(:,k);
                P_k = prob_bc.P(:,:,k);
                constraint_val = alpha' * mu_k + sqrt(alpha' * P_k * alpha) * inv_norm_p;
                violation = constraint_val - beta;
                max_viol_bc = max(max_viol_bc, violation);
            end
        end
    end
    max_viol_bc_str = sprintf('%.6e', max_viol_bc);
end

% Compute FullCovarianceSteering metrics
if ~isempty(prob_fc)
    time_fc_str = sprintf('%.4f', time_fc);
    iters_fc_str = sprintf('%d', length(wrapper_fc.iter_history) - 1);
    J_fc = prob_fc.optimal_objective;
    J_fc_str = sprintf('%.6f', J_fc);
    
    P_term_fc = prob_fc.P(:,:,end);
    % Compute error as largest eigenvalue of (P_term_fc - SigmaN)
    % If largest eigenvalue <= 0, achieved covariance is smaller/equal to target, error = 0
    diff_cov_fc = P_term_fc - SigmaN;
    diff_cov_fc = (diff_cov_fc + diff_cov_fc') / 2; % Ensure symmetry
    max_eig_fc = max(eig(diff_cov_fc));
    err_term_fc = max(0, max_eig_fc);  % Error is 0 if max_eig <= 0
    err_term_fc_str = sprintf('%.6e', err_term_fc);
    
    mu_term_fc = prob_fc.mu(:,end);
    err_mu_fc = norm(mu_term_fc - muN);
    err_mu_fc_str = sprintf('%.6e', err_mu_fc);
    
    % Chance constraint satisfaction
    if ~exist('inv_norm_p', 'var')
        inv_norm_p = norminv(1 - p_violation);
    end
    max_viol_fc = 0;
    for j = 1:length(state_cc)
        cc = state_cc{j};
        if strcmp(cc.type, 'affine')
            alpha = cc.alpha;
            beta = cc.beta;
            for k = 1:N+1
                mu_k = prob_fc.mu(:,k);
                P_k = prob_fc.P(:,:,k);
                constraint_val = alpha' * mu_k + sqrt(alpha' * P_k * alpha) * inv_norm_p;
                violation = constraint_val - beta;
                max_viol_fc = max(max_viol_fc, violation);
            end
        end
    end
    max_viol_fc_str = sprintf('%.6e', max_viol_fc);
end

% Compute SqrtQRCovarianceSteering metrics
if ~isempty(prob_qr)
    time_qr_str = sprintf('%.4f', time_qr);
    iters_qr_str = sprintf('%d', prob_qr.scp.report.iters);
    
    % Compute LQR objective value
    J_qr = 0;
    for k = 1:prob_qr.N
        J_qr = J_qr + trace(Q * prob_qr.P(:,:,k)) + trace(R * prob_qr.P_u(:,:,k)) ...
            + prob_qr.mu(:,k)' * Q * prob_qr.mu(:,k) + prob_qr.v(:,k)' * R * prob_qr.v(:,k);
    end
    J_qr_str = sprintf('%.6f', J_qr);
    
    P_term_qr = prob_qr.P(:,:,end);
    % Compute error as largest eigenvalue of (P_term_qr - SigmaN)
    % If largest eigenvalue <= 0, achieved covariance is smaller/equal to target, error = 0
    diff_cov_qr = P_term_qr - SigmaN;
    diff_cov_qr = (diff_cov_qr + diff_cov_qr') / 2; % Ensure symmetry
    max_eig_qr = max(eig(diff_cov_qr));
    err_term_qr = max(0, max_eig_qr);  % Error is 0 if max_eig <= 0
    err_term_qr_str = sprintf('%.6e', err_term_qr);
    
    mu_term_qr = prob_qr.mu(:,end);
    err_mu_qr = norm(mu_term_qr - muN);
    err_mu_qr_str = sprintf('%.6e', err_mu_qr);
    
    % Chance constraint satisfaction
    if ~exist('inv_norm_p', 'var')
        inv_norm_p = norminv(1 - p_violation);
    end
    max_viol_qr = 0;
    for j = 1:length(state_cc)
        cc = state_cc{j};
        if strcmp(cc.type, 'affine')
            alpha = cc.alpha;
            beta = cc.beta;
            for k = 1:N+1
                mu_k = prob_qr.mu(:,k);
                P_k = prob_qr.P(:,:,k);
                constraint_val = alpha' * mu_k + sqrt(alpha' * P_k * alpha) * inv_norm_p;
                violation = constraint_val - beta;
                max_viol_qr = max(max_viol_qr, violation);
            end
        end
    end
    max_viol_qr_str = sprintf('%.6e', max_viol_qr);
end

% Print all metrics (each metric appears once)
fprintf('%-40s %-25s %-25s %-25s\n', 'Solution time (s)', time_bc_str, time_fc_str, time_qr_str);
fprintf('%-40s %-25s %-25s %-25s\n', 'Number of SCP iterations', 'N/A', iters_fc_str, iters_qr_str);
fprintf('%-40s %-25s %-25s %-25s\n', 'Objective value', J_bc_str, J_fc_str, J_qr_str);
fprintf('%-40s %-25s %-25s %-25s\n', 'Terminal cov. violation (max eig)', err_term_bc_str, err_term_fc_str, err_term_qr_str);
fprintf('%-40s %-25s %-25s %-25s\n', 'Terminal mean violation (2-norm)', err_mu_bc_str, err_mu_fc_str, err_mu_qr_str);
fprintf('%-40s %-25s %-25s %-25s\n', 'Max chance constraint violation', max_viol_bc_str, max_viol_fc_str, max_viol_qr_str);

fprintf('%s\n', repmat('-', 1, 115));
