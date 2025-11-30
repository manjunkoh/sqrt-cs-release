%% Scalability vs Horizon Size: 2D Double Integrator with Cone Constraints
% Based on example_cone_path.m
% Fixed total horizon time; dt changes with N
clear; clc;
addpath ./SCvxStar/src/
addpath(genpath('./utils'))
addpath('./src')

figure_settings

% Fixed total time and horizon sizes to test
T_total = 4.0;                  % total time [s] (N=20 with dt=0.2 gives T=4.0)
N_list = [40];      % different horizon lengths to test
num_trials = 1;                 % trials per N for averaging

% Problem dimensions for 2D double integrator
nx = 4; % [x, y, vx, vy]
nu = 2; % [ax, ay]
nw = 4; % process noise dimension same as state

% SCP parameters for the square root method
scp_params = SCPParams();
scp_params.tol_opt = 1E-3;
scp_params.tol_feas = 1E-4;
scp_params.k_max = 100;

% SDP settings
sdp_settings = sdpsettings('verbose', 0, 'solver', 'mosek');

% Cone-shaped chance constraints: -0.2*(x-1) <= y <= 0.2*(x-1)
% Rearranged to alpha' * x <= beta form for state = [x, y, vx, vy]':
% Constraint 1: y <= 0.2(x-1)   ->  y <= 0.2x - 0.2  ->  0.2x + y <= 0.2
% Constraint 2: y >= -0.2(x-1)  ->  y >= -0.2x + 0.2  ->  0.2x - y <= 0.2
alpha1 = [0.2; 1; 0; 0]; beta1 = 0.2;  % Upper bound: y <= 0.2(x-1)
alpha2 = [0.2; -1; 0; 0]; beta2 = 0.2; % Lower bound: y >= -0.2(x-1)

% Violation probability for chance constraints
p_violation = 0.005;

% Results storage
results = struct();
results.N_list = N_list;
results.sqrtqr_time = zeros(length(N_list), num_trials);
results.fullcov_time = zeros(length(N_list), num_trials);
results.blockcholesky_time = zeros(length(N_list), num_trials);
results.sqrtqr_status = cell(length(N_list), num_trials);
results.fullcov_status = cell(length(N_list), num_trials);
results.blockcholesky_status = cell(length(N_list), num_trials);
results.sqrtqr_opt = NaN(length(N_list), num_trials);
results.fullcov_opt = NaN(length(N_list), num_trials);
results.blockcholesky_opt = NaN(length(N_list), num_trials);
results.sqrtqr_success = false(length(N_list), num_trials);
results.sqrtqr_iters = zeros(length(N_list), num_trials);
results.fullcov_iters = zeros(length(N_list), num_trials);

fprintf('Scalability across horizon sizes with fixed total time T=%.2f s\n', T_total);
fprintf('Horizon sizes: %s\n', mat2str(N_list));
fprintf('Trials per setting: %d\n\n', num_trials);

rng(1)

for iN = 1:length(N_list)
    N = N_list(iN);
    dt = T_total / N;
    fprintf('============================================================\n');
    fprintf('Testing N = %d (dt = %.4f s)\n', N, dt);

    % Build discrete-time double integrator for 2D
    A = [1  0   dt  0;
         0  1   0   dt;
         0  0   1   0;
         0  0   0   1];

    B = [0.5*dt^2  0;
         0         0.5*dt^2;
         dt        0;
         0         dt];

    % Time-invariant replication along horizon
    A_sys = repmat(A, [1, 1, N]);
    B_sys = repmat(B, [1, 1, N]);

    % Process noise (diagonal)
    G = diag([0.01, 0.01, 0.01, 0.01]);
    G_sys = repmat(G, [1, 1, N]);

    % Covariance boundary conditions
    Sigma0 = diag([0.1, 0.1, 0.01, 0.01]);
    SigmaN = 0.5 * Sigma0;

    % Mean boundary conditions
    mu0 = [-10; 1; 0; 0];
    muN = zeros(nx,1);

    % Objective weights
    Q = diag([1e-2, 1e-2, 1e-3, 1e-3]);
    R = diag([1, 1]);

    % Build chance constraint structs (apply to all time steps k=1:N+1)
    state_cc = {
        struct('type', 'affine', 'alpha', alpha1, 'beta', beta1, 'p', p_violation, 'nodes', 1:N+1);
        struct('type', 'affine', 'alpha', alpha2, 'beta', beta2, 'p', p_violation, 'nodes', 1:N+1);
    };

    for trial = 1:num_trials
        fprintf('  Trial %2d/%2d ... ', trial, num_trials);

        % BlockCholeskySteering method
        prob_block = BlockCholeskySteering(...
            A=A_sys, B=B_sys, G=G_sys, ...
            P_0=Sigma0, P_f=SigmaN, ...
            Q=Q, R=R, ...
            N=N, ...
            chance_constraints_state=state_cc, ...
            mu_0=mu0, mu_f=muN);

        t_block = tic;
        diag_block = prob_block.solve(sdp_settings);
        results.blockcholesky_time(iN, trial) = toc(t_block);

        % Interpret status and optimal value for block Cholesky
        switch diag_block.problem
            case 0
                results.blockcholesky_status{iN, trial} = 'SOLVED';
                results.blockcholesky_opt(iN, trial) = prob_block.optimal_objective;
            case 1
                results.blockcholesky_status{iN, trial} = 'INFEASIBLE';
            case 4
                results.blockcholesky_status{iN, trial} = 'NUMERICAL';
                results.blockcholesky_opt(iN, trial) = prob_block.optimal_objective;
            otherwise
                results.blockcholesky_status{iN, trial} = sprintf('ERR%d', diag_block.problem);
        end

        % FullCovarianceSteering (iterative wrapper)
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
            verbose=false);

        t_full = tic;
        [diag_full, prob_full] = wrapper_fc.solve(sdp_settings);
        results.fullcov_time(iN, trial) = toc(t_full);
        results.fullcov_iters(iN, trial) = length(wrapper_fc.iter_history) - 1;

        % Interpret status and optimal value
        switch diag_full.problem
            case 0
                results.fullcov_status{iN, trial} = 'SOLVED';
                results.fullcov_opt(iN, trial) = prob_full.optimal_objective;
            case 1
                results.fullcov_status{iN, trial} = 'INFEASIBLE';
            case 4
                results.fullcov_status{iN, trial} = 'NUMERICAL';
                results.fullcov_opt(iN, trial) = prob_full.optimal_objective;
            otherwise
                results.fullcov_status{iN, trial} = sprintf('ERR%d', diag_full.problem);
        end

        % Initial guess for SqrtQR method
        init = struct();
        init.S = interpolate_lower_triangular(chol(Sigma0, 'lower'), chol(SigmaN, 'lower'), N+1, 'log-cholesky');
        init.L = zeros(nu, nx, N);
        init.mu = linspace_vec(mu0, muN, N+1);
        init.v = zeros(nu, N);

        prob_qr = SqrtQRCovarianceSteeringOptimizer(init, ...
            N=N, ...
            A_sys=A_sys, B_sys=B_sys, G_sys=G_sys, ...
            P_0=Sigma0, P_f=SigmaN, Q=Q, R=R, ...
            objective_type='LQR', ...
            chance_constraints_state=state_cc, ...
            mu_0=mu0, mu_f=muN);

        flag = prob_qr.solve(save_bool=false, scp_params=scp_params, verbose=false);
        results.sqrtqr_time(iN, trial) = seconds(prob_qr.scp.report.time);
        results.sqrtqr_success(iN, trial) = flag;
        results.sqrtqr_status{iN, trial} = SCPCode.toString(prob_qr.scp.report.code);
        results.sqrtqr_iters(iN, trial) = prob_qr.scp.report.iters;
        if flag
            prob_qr.postprocess();
            % Compute LQR objective value
            J_qr = 0;
            for k = 1:prob_qr.N
                J_qr = J_qr + trace(Q * prob_qr.P(:,:,k)) + trace(R * prob_qr.P_u(:,:,k)) ...
                    + prob_qr.mu(:,k)' * Q * prob_qr.mu(:,k) + prob_qr.v(:,k)' * R * prob_qr.v(:,k);
            end
            results.sqrtqr_opt(iN, trial) = J_qr;
        end

        fprintf('block=%s (%.2fs), full=%s (%.2fs, %d iters), qr=%s (%.2fs, %d iters)\n', ...
            results.blockcholesky_status{iN, trial}, results.blockcholesky_time(iN, trial), ...
            results.fullcov_status{iN, trial}, results.fullcov_time(iN, trial), results.fullcov_iters(iN, trial), ...
            results.sqrtqr_status{iN, trial}, results.sqrtqr_time(iN, trial), results.sqrtqr_iters(iN, trial));
    end

    % Summary per N
    valid_qr = results.sqrtqr_time(iN, results.sqrtqr_time(iN, :) > 0);
    if ~isempty(valid_qr)
        fprintf('  QR avg time: %.2f s\n', mean(valid_qr));
    end
    fprintf('  FullCov avg time: %.2f s\n', mean(results.fullcov_time(iN, :)));
    fprintf('  BlockCholesky avg time: %.2f s\n\n', mean(results.blockcholesky_time(iN, :)));
end

% Console table summary
fprintf('\n=== TERMINATION AND TIME SUMMARY (Fixed T=%.2fs) ===\n', T_total);
fprintf(' N  |   dt    | BlockChol | FullCov | SqrtQR  | BlockCholOptVal | FullCovOptVal | SqrtQROptVal | BlockTime(s) | FullTime(s) | QRTime(s)\n');
fprintf('----|---------|-----------|---------|---------|-----------------|---------------|--------------|--------------|-------------|----------\n');
for iN = 1:length(N_list)
    N = N_list(iN);
    dt = T_total / N;
    for trial = 1:num_trials
        block_opt_str = 'N/A';
        if ~isnan(results.blockcholesky_opt(iN, trial))
            block_opt_str = sprintf('%.2e', results.blockcholesky_opt(iN, trial));
        end
        cov_opt_str = 'N/A';
        if ~isnan(results.fullcov_opt(iN, trial))
            cov_opt_str = sprintf('%.2e', results.fullcov_opt(iN, trial));
        end
        qr_opt_str = 'N/A';
        if ~isnan(results.sqrtqr_opt(iN, trial))
            qr_opt_str = sprintf('%.2e', results.sqrtqr_opt(iN, trial));
        end
        fprintf('%3d | %7.4f | %-9s | %-7s | %-7s | %15s | %13s | %12s | %12.2f | %11.2f | %8.2f\n', ...
            N, dt, ...
            results.blockcholesky_status{iN, trial}, ...
            results.fullcov_status{iN, trial}, ...
            results.sqrtqr_status{iN, trial}, ...
            block_opt_str, cov_opt_str, qr_opt_str, ...
            results.blockcholesky_time(iN, trial), ...
            results.fullcov_time(iN, trial), ...
            results.sqrtqr_time(iN, trial));
    end
    if iN < length(N_list)
        fprintf('----|---------|-----------|---------|---------|-----------------|---------------|--------------|--------------|-------------|----------\n');
    end
end

%% Simple plot: average times vs N
avg_full = mean(results.fullcov_time, 2);
avg_block = mean(results.blockcholesky_time, 2);
avg_qr = zeros(length(N_list), 1);
for iN = 1:length(N_list)
    times = results.sqrtqr_time(iN, results.sqrtqr_time(iN, :) > 0);
    if ~isempty(times)
        avg_qr(iN) = mean(times);
    else
        avg_qr(iN) = NaN;
    end
end

% Compute average objective values
avg_full_opt = zeros(length(N_list), 1);
avg_block_opt = zeros(length(N_list), 1);
avg_qr_opt = zeros(length(N_list), 1);
for iN = 1:length(N_list)
    % Full covariance: average over non-NaN values
    valid_full = results.fullcov_opt(iN, ~isnan(results.fullcov_opt(iN, :)));
    if ~isempty(valid_full)
        avg_full_opt(iN) = mean(valid_full);
    else
        avg_full_opt(iN) = NaN;
    end
    
    % Block Cholesky: average over non-NaN values
    valid_block = results.blockcholesky_opt(iN, ~isnan(results.blockcholesky_opt(iN, :)));
    if ~isempty(valid_block)
        avg_block_opt(iN) = mean(valid_block);
    else
        avg_block_opt(iN) = NaN;
    end
    
    % SqrtQR: average over non-NaN values
    valid_qr = results.sqrtqr_opt(iN, ~isnan(results.sqrtqr_opt(iN, :)));
    if ~isempty(valid_qr)
        avg_qr_opt(iN) = mean(valid_qr);
    else
        avg_qr_opt(iN) = NaN;
    end
end

%% Plot runtime
figure(Position=[0, 0, 20, 12]);
loglog(N_list, avg_block, '^-', 'Color', '#D55E00', 'LineWidth', 2, 'MarkerSize', 8); hold on;
loglog(N_list, avg_full, 'o-', 'Color', '#0082B2', 'LineWidth', 2, 'MarkerSize', 8); hold on;
loglog(N_list, avg_qr, 's-', 'Color', '#000000', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Horizon length N');
ylabel('Average runtime (s)');
xticks(N_list)
legend('Okamoto \& Tsiotras (2019)', 'Liu et al. (2025)', 'Proposed method', 'Location', 'northwest', 'EdgeColor', 'none');
exportgraphics(gcf, 'figures/horizon_size_scalability_constrained.png', Resolution=300)
exportgraphics(gcf, 'figures/horizon_size_scalability_constrained.pdf', ContentType='vector')

%% Plot objective function values
figure(Position=[0, 0, 20, 12]);
semilogy(N_list, avg_block_opt./avg_full_opt, '^-', 'Color', '#D55E00', 'LineWidth', 2, 'MarkerSize', 8); hold on;
semilogy(N_list, avg_qr_opt./avg_full_opt, 's-', 'Color', '#000000', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Horizon length N');
ylabel('Cost ratio to Liu et al.');
xticks(N_list)
lgd = legend('Okamoto \& Tsiotras (2019)', 'Proposed method', 'Location', 'west', 'EdgeColor', 'none');
xlim([N_list(1) N_list(end)])
if any(~isnan(avg_block_opt./avg_full_opt))
    ylim([0.99, max(avg_block_opt./avg_full_opt)])
end

exportgraphics(gcf, 'figures/horizon_size_cost_comparison_constrained.png', Resolution=300)
exportgraphics(gcf, 'figures/horizon_size_cost_comparison_constrained.pdf', ContentType='vector')

%%
save('data/horizon_size_scalability_constrained_results.mat', 'results', 'T_total', 'N_list');
fprintf('\nSaved results to data/horizon_size_scalability_constrained_results.mat\n');

