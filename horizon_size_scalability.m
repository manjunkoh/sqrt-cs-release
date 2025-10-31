%% Scalability vs Horizon Size: 3D Double Integrator System
% Fixed total horizon time; dt changes with N
clear;
clc;
addpath ../SCvxStar/src/
addpath(genpath('./utils'))

figure_settings

% Fixed total time and horizon sizes to test
T_total = 3.0;                  % total time [s]
N_list = [10, 20, 40, 80];      % different horizon lengths to test
num_trials = 1;                % trials per N for averaging

% Problem dimensions for 3D double integrator
nx = 6; % [px, py, pz, vx, vy, vz]
nu = 3; % [ax, ay, az]
nw = nx;

% SCP parameters for the square root method
scp_params = SCPParams();
scp_params.tol_opt = 1E-2;
scp_params.tol_feas = 1E-4;
scp_params.k_max = 200;

% Results storage
results = struct();
results.N_list = N_list;
results.sqrtqr_time = zeros(length(N_list), num_trials);
results.fullcov_time = zeros(length(N_list), num_trials);
results.sqrtqr_status = cell(length(N_list), num_trials);
results.fullcov_status = cell(length(N_list), num_trials);
results.sqrtqr_opt = NaN(length(N_list), num_trials);
results.fullcov_opt = NaN(length(N_list), num_trials);
results.sqrtqr_success = false(length(N_list), num_trials);

fprintf('Scalability across horizon sizes with fixed total time T=%.2f s\n', T_total);
fprintf('Horizon sizes: %s\n', mat2str(N_list));
fprintf('Trials per setting: %d\n\n', num_trials);

rng(1)

for iN = 1:length(N_list)
    N = N_list(iN);
    dt = T_total / N;
    fprintf('============================================================\n');
    fprintf('Testing N = %d (dt = %.4f s)\n', N, dt);

    % Build discrete-time double integrator for 3D
    I3 = eye(3);
    Z3 = zeros(3);
    A = [I3, dt * I3; Z3, I3];
    B = [0.5 * dt^2 * I3; dt * I3];

    % Time-invariant replication along horizon
    A_sys = repmat(A, [1, 1, N]);
    B_sys = repmat(B, [1, 1, N]);

    % Process noise matrix (scaled with dt for stability of covariance growth)
    q = 0.05; % base spectral density
    G = sqrt(q * dt) * eye(nw);
    G_sys = repmat(G, [1, 1, N]);

    % Covariance boundary conditions
    Sigma0 = eye(nx);
    SigmaN = 0.5 * Sigma0;

    % Objective weights
    Q = eye(nx);
    R = 10 * eye(nu);

    for trial = 1:num_trials
        fprintf('  Trial %2d/%2d ... ', trial, num_trials);

        % Full covariance method
        prob_full = CovarianceSteering( ...
            A=A_sys, B=B_sys, D=G_sys, ...
            P_0=Sigma0, P_f=SigmaN, ...
            Q=Q, R=R, ...
            nx=nx, nu=nu, N=N);

        t_full = tic;
        diag_full = prob_full.solve_problem(verbose=2);
        results.fullcov_time(iN, trial) = toc(t_full);

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

        % If clearly problematic (not numerical), skip QR to keep stats clean
        if diag_full.problem && diag_full.problem ~= 4
            fprintf('full=%s (%.2fs), skip QR\n', results.fullcov_status{iN, trial}, results.fullcov_time(iN, trial));
            continue;
        end

        % Initial guess for SqrtQR method
        init.S = interpolate_lower_triangular(chol(Sigma0, 'lower'), chol(SigmaN, 'lower'), N+1, 'log-cholesky');
        init.L = zeros(nu, nx, N);
        init.mu = zeros(nx, N+1);
        init.v = zeros(nu, N);

        prob = SqrtQRCovarianceSteering(init, ...
            N=N, nx=nx, nu=nu, nw=nw, ...
            A_sys=A_sys, B_sys=B_sys, G_sys=G_sys, ...
            P_0=Sigma0, P_f=SigmaN, Q=Q, R=R, ...
            objective_type='LQR');

        flag = prob.solve(save_bool=false, scp_params=scp_params, verbose=true);
        results.sqrtqr_time(iN, trial) = seconds(prob.scp.report.time);
        results.sqrtqr_success(iN, trial) = flag;
        results.sqrtqr_status{iN, trial} = SCPCode.toString(prob.scp.report.code);
        if flag
            results.sqrtqr_opt(iN, trial) = prob.optimal_objective;
        end

        fprintf('full=%s (%.2fs), qr=%s (%.2fs)\n', ...
            results.fullcov_status{iN, trial}, results.fullcov_time(iN, trial), ...
            results.sqrtqr_status{iN, trial}, results.sqrtqr_time(iN, trial));
    end

    % Summary per N
    valid_qr = results.sqrtqr_time(iN, results.sqrtqr_time(iN, :) > 0);
    if ~isempty(valid_qr)
        fprintf('  QR avg time: %.2f s\n', mean(valid_qr));
    end
    fprintf('  FullCov avg time: %.2f s\n\n', mean(results.fullcov_time(iN, :)) );
end

% Console table summary
fprintf('\n=== TERMINATION AND TIME SUMMARY (Fixed T=%.2fs) ===\n', T_total);
fprintf(' N  |   dt    | FullCov Status | SqrtQR Status  | FullCovOptVal | SqrtQROptVal | FullTime(s) | QRTime(s)\n');
fprintf('----|---------|----------------|----------------|---------------|--------------|-------------|----------\n');
for iN = 1:length(N_list)
    N = N_list(iN);
    dt = T_total / N;
    for trial = 1:num_trials
        cov_opt_str = 'N/A';
        if ~isnan(results.fullcov_opt(iN, trial))
            cov_opt_str = sprintf('%.2e', results.fullcov_opt(iN, trial));
        end
        qr_opt_str = 'N/A';
        if ~isnan(results.sqrtqr_opt(iN, trial))
            qr_opt_str = sprintf('%.2e', results.sqrtqr_opt(iN, trial));
        end
        fprintf('%3d | %7.4f | %-14s | %-14s | %13s | %12s | %10.2f | %8.2f\n', ...
            N, dt, ...
            results.fullcov_status{iN, trial}, ...
            results.sqrtqr_status{iN, trial}, ...
            cov_opt_str, qr_opt_str, ...
            results.fullcov_time(iN, trial), ...
            results.sqrtqr_time(iN, trial));
    end
    if iN < length(N_list)
        fprintf('----|---------|----------------|----------------|---------------|--------------|------------|----------\n');
    end
end

%% Simple plot: average times vs N
avg_full = mean(results.fullcov_time, 2);
avg_qr = zeros(length(N_list), 1);
for iN = 1:length(N_list)
    times = results.sqrtqr_time(iN, results.sqrtqr_time(iN, :) > 0);
    if ~isempty(times)
        avg_qr(iN) = mean(times);
    else
        avg_qr(iN) = NaN;
    end
end

figure();
loglog(N_list, avg_full, 'o-b', 'LineWidth', 2, 'MarkerSize', 8); hold on;
loglog(N_list, avg_qr, 's-r', 'LineWidth', 2, 'MarkerSize', 8);
grid on;
xlabel('Horizon length N');
ylabel('Average runtime (s)');
title(sprintf('Scalability vs N (T=%.2fs, 3D double integrator)', T_total));
legend('Full Covariance', 'SqrtQR', 'Location', 'northwest');

save('data/horizon_size_scalability_results.mat', 'results', 'T_total', 'N_list');
fprintf('\nSaved results to /data/horizon_size_scalability_results.mat\n');

