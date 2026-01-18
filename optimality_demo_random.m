%% Numerical demonstration of global optimality of sqrt method with Random System Dynamics
% Uses MATLAB's drss command to generate random discrete-time systems

clear; clc; clear all
addpath ./SCvxStar/src/
addpath(genpath('./utils'))
addpath('./src')

figure_settings

%%
% Fixed horizon length
N = 30; % horizon length
num_trials = 20;                % number of trials for averaging

% Problem dimensions for random system
nx = 6;  % State dimension
nu = 3;  % Input dimension (nx/2)
nw = nx; % Noise dimension

% SCP parameters for the square root method
scp_params = SCPParams();
scp_params.tol_opt = 1E-5;
scp_params.tol_feas = 1E-5;
scp_params.k_max = 300;
% scp_params.r_init = 1.0;
scp_params.r_init = 2.0;
% scp_params.superlinear = true;

% Results storage
results = struct();
results.N = N;
results.sqrtqr_time = zeros(1, num_trials);
results.fullcov_time = zeros(1, num_trials);
results.fullcov_solvertime = zeros(1, num_trials);
results.fullcov_yalmiptime = zeros(1, num_trials);
results.blockcholesky_time = zeros(1, num_trials);
results.blockcholesky_solvertime = zeros(1, num_trials);
results.blockcholesky_yalmiptime = zeros(1, num_trials);
results.sqrtqr_status = cell(1, num_trials);
results.fullcov_status = cell(1, num_trials);
results.blockcholesky_status = cell(1, num_trials);
results.sqrtqr_opt = NaN(1, num_trials);
results.fullcov_opt = NaN(1, num_trials);
results.blockcholesky_opt = NaN(1, num_trials);
results.sqrtqr_success = false(1, num_trials);
% Storage for system matrices and problem objects
results.A = cell(1, num_trials);
results.B = cell(1, num_trials);
results.G = cell(1, num_trials);
results.prob = cell(1, num_trials);
results.prob_full = cell(1, num_trials);

fprintf('State size: nx=%d, nu=%d, nw=%d\n', nx, nu, nw);
fprintf('Horizon size: %d\n', N);
fprintf('Trials %d\n\n', num_trials);

%%
rng(1)

fprintf('============================================================\n');
fprintf('Testing N = %d\n', N);

for trial = 1:num_trials
        fprintf('  Trial %2d/%2d ... ', trial, num_trials);
        
        % Generate random discrete-time system using drss
        sys = drss(nx, nx, nu);
        A = sys.A;
        B = sys.B;

        % Time-invariant replication along horizon
        A_sys = repmat(A, [1, 1, N]);
        B_sys = repmat(B, [1, 1, N]);
        assert(rank(ctrb(sys)) == nx)

        % Process noise matrix
        q = 0.05; % base spectral density
        G = sqrt(q) * eye(nw);
        G_sys = repmat(G, [1, 1, N]);
        
        % Store system matrices for this trial
        results.A{trial} = A;
        results.B{trial} = B;
        results.G{trial} = G;

        % Covariance boundary conditions
        Sigma0 = eye(nx);
        SigmaN = 0.5 * Sigma0;

        % Objective weights
        Q = 0.1 * eye(nx);
        R = eye(nu);

        % Full covariance method
        prob_full = FullCovarianceSteering( ...
            A=A_sys, B=B_sys, G=G_sys, ...
            P_0=Sigma0, P_f=SigmaN, ...
            Q=Q, R=R, ...
            N=N);
        
        % Store prob_full for this trial
        results.prob_full{trial} = prob_full;

        t_full = tic;
        diag_full = prob_full.solve();
        results.fullcov_time(trial) = toc(t_full);
        
        % Store solver and YALMIP times from diagnostics
        if isfield(diag_full, 'solvertime')
            results.fullcov_solvertime(trial) = diag_full.solvertime;
        end
        if isfield(diag_full, 'yalmiptime')
            results.fullcov_yalmiptime(trial) = diag_full.yalmiptime;
        end

        % Interpret status and optimal value
        switch diag_full.problem
            case 0
                results.fullcov_status{trial} = 'SOLVED';
                results.fullcov_opt(trial) = prob_full.optimal_objective;
            case 1
                results.fullcov_status{trial} = 'INFEASIBLE';
            case 4
                results.fullcov_status{trial} = 'NUMERICAL';
                results.fullcov_opt(trial) = prob_full.optimal_objective;
            otherwise
                results.fullcov_status{trial} = sprintf('ERR%d', diag_full.problem);
        end

        % Block Cholesky method
        % prob_block = BlockCholeskySteering( ...
        %     A=A_sys, B=B_sys, G=G_sys, ...
        %     P_0=Sigma0, P_f=SigmaN, ...
        %     Q=Q, R=R, ...
        %     N=N);
        % 
        % t_block = tic;
        % diag_block = prob_block.solve();
        % results.blockcholesky_time(iN, trial) = toc(t_block);
        % 
        % % Store solver and YALMIP times from diagnostics
        % if isfield(diag_block, 'solvertime')
        %     results.blockcholesky_solvertime(iN, trial) = diag_block.solvertime;
        % end
        % if isfield(diag_block, 'yalmiptime')
        %     results.blockcholesky_yalmiptime(iN, trial) = diag_block.yalmiptime;
        % end
        % 
        % % Interpret status and optimal value for block Cholesky
        % switch diag_block.problem
        %     case 0
        %         results.blockcholesky_status{iN, trial} = 'SOLVED';
        %         results.blockcholesky_opt(iN, trial) = prob_block.optimal_objective;
        %     case 1
        %         results.blockcholesky_status{iN, trial} = 'INFEASIBLE';
        %     case 4
        %         results.blockcholesky_status{iN, trial} = 'NUMERICAL';
        %         results.blockcholesky_opt(iN, trial) = prob_block.optimal_objective;
        %     otherwise
        %         results.blockcholesky_status{iN, trial} = sprintf('ERR%d', diag_block.problem);
        % end

        % If clearly problematic (not numerical), skip QR to keep stats clean
        % if diag_full.problem && diag_full.problem ~= 4
        %     fprintf('full=%s (%.2fs), block=%s (%.2fs), skip QR\n', ...
        %         results.fullcov_status{iN, trial}, results.fullcov_time(iN, trial), ...
        %         results.blockcholesky_status{iN, trial}, results.blockcholesky_time(iN, trial));
        %     continue;
        % end

        % Initial guess for SqrtQR method
        init = struct();
        init.S = interpolate_lower_triangular(chol(Sigma0, 'lower'), chol(SigmaN, 'lower'), N+1, 'log-cholesky');
        init.L = zeros(nu, nx, N);
        init.mu = zeros(nx, N+1);
        init.v = zeros(nu, N);

        prob = SqrtQRCovarianceSteering(init, ...
            N=N, ...
            A_sys=A_sys, B_sys=B_sys, G_sys=G_sys, ...
            P_0=Sigma0, P_f=SigmaN, Q=Q, R=R, ...
            objective_type='LQR');

        flag = prob.solve(save_bool=false, scp_params=scp_params, verbose=true);
        if flag
            prob.postprocess();
        end
        results.sqrtqr_time(trial) = seconds(prob.scp.report.time);
        results.sqrtqr_success(trial) = flag;
        results.sqrtqr_status{trial} = SCPCode.toString(prob.scp.report.code);
        if flag
            results.sqrtqr_opt(trial) = prob.optimal_objective;
        end
        
        % Store prob for this trial
        results.prob{trial} = prob;

        fprintf('full=%s (%.2fs), block=%s (%.2fs), qr=%s (%.2fs)\n', ...
            results.fullcov_status{trial}, results.fullcov_time(trial), ...
            results.blockcholesky_status{trial}, results.blockcholesky_time(trial), ...
            results.sqrtqr_status{trial}, results.sqrtqr_time(trial));
end

% Summary
valid_qr = results.sqrtqr_time(results.sqrtqr_time > 0);
if ~isempty(valid_qr)
    fprintf('  QR avg time: %.2f s\n', mean(valid_qr));
end
fprintf('  FullCov avg time: %.2f s\n', mean(results.fullcov_time));
% fprintf('  BlockCholesky avg time: %.2f s\n\n', mean(results.blockcholesky_time));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Console table summary
fprintf('\n=== TERMINATION AND TIME SUMMARY (Random Dynamics) ===\n');
fprintf('Trial | FullCov | BlockChol | SqrtQR  | FullCovOptVal | BlockCholOptVal | SqrtQROptVal | FullTime(s) | BlockTime(s) | QRTime(s)\n');
fprintf('------|---------|-----------|---------|---------------|-----------------|--------------|-------------|--------------|----------\n');
for trial = 1:num_trials
    cov_opt_str = 'N/A';
    if ~isnan(results.fullcov_opt(trial))
        cov_opt_str = sprintf('%.2e', results.fullcov_opt(trial));
    end
    block_opt_str = 'N/A';
    if ~isnan(results.blockcholesky_opt(trial))
        block_opt_str = sprintf('%.2e', results.blockcholesky_opt(trial));
    end
    qr_opt_str = 'N/A';
    if ~isnan(results.sqrtqr_opt(trial))
        qr_opt_str = sprintf('%.2e', results.sqrtqr_opt(trial));
    end
    fprintf('%5d | %-7s | %-9s | %-7s | %13s | %15s | %12s | %10.2f | %12.2f | %8.2f\n', ...
        trial, ...
        results.fullcov_status{trial}, ...
        results.blockcholesky_status{trial}, ...
        results.sqrtqr_status{trial}, ...
        cov_opt_str, block_opt_str, qr_opt_str, ...
        results.fullcov_time(trial), ...
        results.blockcholesky_time(trial), ...
        results.sqrtqr_time(trial));
end

%% Compute statistics
% Average times
avg_full = mean(results.fullcov_time);
avg_block = mean(results.blockcholesky_time);
min_full = min(results.fullcov_time);
max_full = max(results.fullcov_time);
min_block = min(results.blockcholesky_time);
max_block = max(results.blockcholesky_time);

% SqrtQR: average, min and max over valid trials
times = results.sqrtqr_time(results.sqrtqr_time > 0);
if ~isempty(times)
    avg_qr = mean(times);
    min_qr = min(times);
    max_qr = max(times);
else
    avg_qr = NaN;
    min_qr = NaN;
    max_qr = NaN;
end

% Compute average objective values
% Full covariance: average over non-NaN values
valid_full = results.fullcov_opt(~isnan(results.fullcov_opt));
if ~isempty(valid_full)
    avg_full_opt = mean(valid_full);
else
    avg_full_opt = NaN;
end

% Block Cholesky: average over non-NaN values
valid_block = results.blockcholesky_opt(~isnan(results.blockcholesky_opt));
if ~isempty(valid_block)
    avg_block_opt = mean(valid_block);
else
    avg_block_opt = NaN;
end

% SqrtQR: average over non-NaN values
valid_qr = results.sqrtqr_opt(~isnan(results.sqrtqr_opt));
if ~isempty(valid_qr)
    avg_qr_opt = mean(valid_qr);
else
    avg_qr_opt = NaN;
end

%% Plot runtime comparison
figure(Position=[0, 0, 11, 12]);
hold on;

% Create bar plot with error bars
methods = {'Liu et al.', 'Okamoto \& Tsiotras', 'Proposed method'};
avg_times = [avg_full, avg_block, avg_qr];
min_times = [min_full, min_block, min_qr];
max_times = [max_full, max_block, max_qr];

% Plot bars
b = bar(avg_times, 'FaceColor', 'flat');
b.CData(1,:) = [0, 130/255, 178/255]; % #0082B2
b.CData(2,:) = [213/255, 94/255, 0];  % #D55E00
b.CData(3,:) = [0, 0, 0];              % #000000

% Add error bars (min-max range)
for i = 1:length(avg_times)
    if ~isnan(avg_times(i))
        errorbar(i, avg_times(i), avg_times(i) - min_times(i), max_times(i) - avg_times(i), ...
            'k', 'LineWidth', 1.5, 'CapSize', 10);
    end
end

set(gca, 'XTickLabel', methods);
set(gca, 'YScale', 'log');
grid on;
ylabel('Runtime (s)');
title(sprintf('Runtime Comparison (N=%d, Random System Dynamics)', N));
legend('off');
exportgraphics(gcf, 'figures/horizon_size_scalability_random_small.png', Resolution=300)
exportgraphics(gcf, 'figures/horizon_size_scalability_random_small.pdf', ContentType='vector')

%% Plot objective function values
figure(Position=[0, 0, 11, 12]);
methods_cost = {'Okamoto \& Tsiotras', 'Proposed method'};
cost_ratios = [avg_block_opt./avg_full_opt, avg_qr_opt./avg_full_opt];

b = bar(cost_ratios, 'FaceColor', 'flat');
b.CData(1,:) = [213/255, 94/255, 0];  % #D55E00
b.CData(2,:) = [0, 0, 0];              % #000000

set(gca, 'XTickLabel', methods_cost);
grid on;
ylabel('Cost ratio to Liu et al.');
title(sprintf('Cost Comparison (N=%d, Random System Dynamics)', N));
if any(~isnan(cost_ratios))
    ylim([0.99, max(cost_ratios)])
end

exportgraphics(gcf, 'figures/horizon_size_cost_comparision_random_small.png', Resolution=300)
exportgraphics(gcf, 'figures/horizon_size_cost_comparision_random_small.pdf', ContentType='vector')

%%
save('data/horizon_size_scalability_random_results.mat', 'results', 'N', 'nx', 'nu', 'nw');
fprintf('\nSaved results to data/horizon_size_scalability_random_results.mat\n');

