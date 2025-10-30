%% Runtime Comparison: Square Root Covariance Steering with Different State Sizes
% Uses MATLAB's drss command to generate random discrete-time systems
% Compares nx = 4, 8, 16, 32 with nu = nx/2, nw = nx
clear; 
clc;
addpath ../SCvxStar/src/
addpath(genpath('./utils'))

figure_settings

% Problem parameters
N = 30;  % Fixed time horizon
state_sizes = [4];  % Different state dimensions to test
num_trials = 20;  % Number of trials per state size for averaging

% SCP parameters for the square root method
scp_params = SCPParams();
scp_params.tol_opt = 1E-2;
scp_params.tol_feas = 1E-4;
scp_params.k_max = 200;  % Limit iterations for timing

% Storage for results
results = struct(); 
results.state_sizes = state_sizes;
results.runtimes = zeros(length(state_sizes), num_trials);
results.success_flags = false(length(state_sizes), num_trials);
results.avg_runtimes = zeros(length(state_sizes), 1);
results.std_runtimes = zeros(length(state_sizes), 1);

% Storage for termination status tracking
results.covariance_steering_status = cell(length(state_sizes), num_trials);
results.sqrt_qr_status = cell(length(state_sizes), num_trials);
results.covariance_steering_codes = zeros(length(state_sizes), num_trials);
results.sqrt_qr_codes = zeros(length(state_sizes), num_trials);

% Storage for optimal values
results.covariance_steering_optimal = zeros(length(state_sizes), num_trials);
results.sqrt_qr_optimal = zeros(length(state_sizes), num_trials);

fprintf('Starting runtime comparison across different state sizes...\n');
fprintf('State sizes: %s\n', mat2str(state_sizes));
fprintf('Time horizon: %d\n', N);
fprintf('Number of trials per size: %d\n\n', num_trials);

rng(3)

for i = 1:length(state_sizes)
    nx = state_sizes(i);
    nu = nx / 2;  % Input size is half of state size
    nw = nx;  % Noise size is state size
    
    fprintf('=============================================================\n')
    fprintf('Testing nx = %d (nu = %d, nw = %d)...\n', nx, nu, nw);
    
    for trial = 1:num_trials
        fprintf('---------------------------------------------------------\n')
        fprintf('  Trial %d/%d: ', trial, num_trials);
        
        % Generate random discrete-time system using drss
        % drss(nx, nu, ny) where ny is output dimension (we'll use nx)
        sys = drss(nx, nx, nu);
        A = sys.A;
        B = sys.B;
        
        % Ensure system is stable (eigenvalues inside unit circle)
        % if max(abs(eig(A))) >= 1
            % error('Non stable system.')
        % end
        
        % Replicate for time-varying arrays
        A_sys = repmat(A, [1,1,N]);
        B_sys = repmat(B, [1,1,N]);
        
        % Process noise matrix (diagonal with small values)
        G = 0.1 * eye(nw);
        G_sys = repmat(G, [1,1,N]);
        
        % Covariance boundary conditions
        Sigma0 = eye(nx);
        SigmaN = 0.5 * Sigma0;
        
        % Mean boundary conditions
        mu0 = zeros(nx,1);
        muN = zeros(nx,1);
        
        % Objective weights
        Q = eye(nx);  % Identity state cost
        R = 1e3 * eye(nu);  % Control cost

        % Create Covariance Steering Problem
        prob_full_covariance = CovarianceSteering(...
            A=A_sys, B=B_sys, D=G_sys, ...
            P_0=Sigma0, P_f=SigmaN, ...
            Q=Q, R=R, ...
            nx=nx, nu=nu, N=N);

        diagnostic_full_covariance = prob_full_covariance.solve_problem();
        
        % Store CovarianceSteering status and optimal value
        results.covariance_steering_codes(i, trial) = diagnostic_full_covariance.problem;
        switch diagnostic_full_covariance.problem
            case 0
            results.covariance_steering_status{i, trial} = 'SOLVED';
            results.covariance_steering_optimal(i, trial) = prob_full_covariance.optimal_objective;
            case 1
            results.covariance_steering_status{i, trial} = 'INFEASIBLE';
            results.covariance_steering_optimal{i, trial} = NaN;
            case 4
            results.covariance_steering_status{i, trial} = 'NUMERICAL';
            results.covariance_steering_optimal(i, trial) = prob_full_covariance.optimal_objective;
            otherwise
            results.covariance_steering_status{i, trial} = sprintf('Error %d', diagnostic_full_covariance.problem);
            results.covariance_steering_optimal(i, trial) = NaN;
        end
        
        if diagnostic_full_covariance.problem && diagnostic_full_covariance.problem ~= 4
            error('Full covariance method had issues.');
        end

        % Initial guess
        init.S = linspace_mat(chol(Sigma0, 'lower'), chol(SigmaN, 'lower'), N+1);
        init.L = zeros(nu, nx, N);
        init.mu = zeros(nx, N+1);
        init.v = zeros(nu, N);
        
        % Create problem
        prob = SqrtQRCovarianceSteering(init, ...
            N=N, nx=nx, nu=nu, nw=nw, ...
            A_sys=A_sys, ...
            B_sys=B_sys, ...
            G_sys=G_sys, ...
            P_0=Sigma0, P_f=SigmaN, Q=Q, R=R, ...
            objective_type='LQR');
        
        % Time the solve
        flag_solved = prob.solve(save_bool=false, scp_params=scp_params, verbose=true);
        
        % Store SqrtQRCovarianceSteering status and optimal value
        results.sqrt_qr_codes(i, trial) = prob.scp.report.code;
        results.sqrt_qr_status{i, trial} = SCPCode.toString(prob.scp.report.code);
        
        if flag_solved
            results.sqrt_qr_optimal(i, trial) = prob.optimal_objective;
        else
            results.sqrt_qr_optimal(i, trial) = NaN;
        end
        
        results.runtimes(i, trial) = seconds(prob.scp.report.time);
        results.success_flags(i, trial) = flag_solved;
        
        if flag_solved
            % fprintf('QR Square Root method solved in %.2f seconds\n', solve_time);
        else
            fprintf('QR Square Root method terminated with status %s in %.2f seconds\n', ...
                SCPCode.toString(prob.scp.report.code), seconds(prob.scp.report.time));
        end
            
    end
    
    % Calculate statistics for this state size
    valid_times = results.runtimes(i, ~isnan(results.runtimes(i, :)));
    if ~isempty(valid_times)
        results.avg_runtimes(i) = mean(valid_times);
        results.std_runtimes(i) = std(valid_times);
    else
        results.avg_runtimes(i) = NaN;
        results.std_runtimes(i) = NaN;
    end
    
    fprintf('  Average runtime: %.2f ± %.2f seconds\n\n', ...
        results.avg_runtimes(i), results.std_runtimes(i));
end

% Display summary results
fprintf('\n=== RUNTIME COMPARISON SUMMARY ===\n');
fprintf('State Size | Avg Runtime (s) | Std Dev (s) | Success Rate\n');
fprintf('-----------|-----------------|-------------|-------------\n');
for i = 1:length(state_sizes)
    success_rate = sum(results.success_flags(i, :)) / num_trials * 100;
    fprintf('    %2d     |     %8.2f    |   %8.2f   |    %6.1f%%\n', ...
        state_sizes(i), results.avg_runtimes(i), results.std_runtimes(i), success_rate);
end

% Create status comparison table
fprintf('\n=== TERMINATION STATUS COMPARISON ===\n');
fprintf('Trial | Full Cov Status | SqrtQR Status | FullCovOptVal | SqrtQROptVal | Runtime (s)\n');
fprintf('------|-----------------|---------------|---------------|-------------|-------------\n');

for i = 1:length(state_sizes)
    for trial = 1:num_trials
        % Format optimal values
        cov_opt_str = 'N/A';
        if ~isnan(results.covariance_steering_optimal(i, trial))
            cov_opt_str = sprintf('%.2e', results.covariance_steering_optimal(i, trial));
        end
        
        sqrt_qr_opt_str = 'N/A';
        if ~isnan(results.sqrt_qr_optimal(i, trial))
            sqrt_qr_opt_str = sprintf('%.2e', results.sqrt_qr_optimal(i, trial));
        end
        
        fprintf('  %2d  |      %-13s |    %-11s |    %-12s |   %-9s |    %8.2f\n', ...
            trial, ...
            results.covariance_steering_status{i, trial}, ...
            results.sqrt_qr_status{i, trial}, ...
            cov_opt_str, ...
            sqrt_qr_opt_str, ...
            results.runtimes(i, trial));
    end
    if i < length(state_sizes)
        fprintf('------|-----------------|---------------|---------------|-------------|-------------\n');
    end
end

return
%%
% Create visualization plots
figure('Position', [100, 100, 1200, 800]);

% Subplot 1: Status comparison heatmap
subplot(2,2,1);
status_data = zeros(length(state_sizes), num_trials);
for i = 1:length(state_sizes)
    for trial = 1:num_trials
        % Convert status to numeric for heatmap
        if strcmp(results.covariance_steering_status{i, trial}, 'Success')
            status_data(i, trial) = 1;
        elseif strcmp(results.covariance_steering_status{i, trial}, 'Infeasible')
            status_data(i, trial) = 0.5;
        else
            status_data(i, trial) = 0;
        end
    end
end
imagesc(status_data);
colorbar;
colormap([1 0 0; 1 1 0; 0 1 0]); % Red, Yellow, Green
title('CovarianceSteering Status Heatmap');
xlabel('Trial');
ylabel('State Size');
set(gca, 'XTick', 1:num_trials);
set(gca, 'YTick', 1:length(state_sizes), 'YTickLabel', state_sizes);

% Subplot 2: SqrtQR status heatmap
subplot(2,2,2);
qr_status_data = zeros(length(state_sizes), num_trials);
for i = 1:length(state_sizes)
    for trial = 1:num_trials
        % Convert SCP status to numeric
        if strcmp(results.sqrt_qr_status{i, trial}, 'Solved')
            qr_status_data(i, trial) = 1;
        elseif contains(results.sqrt_qr_status{i, trial}, 'Infeasible')
            qr_status_data(i, trial) = 0.5;
        else
            qr_status_data(i, trial) = 0;
        end
    end
end
imagesc(qr_status_data);
colorbar;
colormap([1 0 0; 1 1 0; 0 1 0]); % Red, Yellow, Green
title('SqrtQRCovarianceSteering Status Heatmap');
xlabel('Trial');
ylabel('State Size');
set(gca, 'XTick', 1:num_trials);
set(gca, 'YTick', 1:length(state_sizes), 'YTickLabel', state_sizes);

% Subplot 3: Runtime comparison
subplot(2,2,3);
errorbar(state_sizes, results.avg_runtimes, results.std_runtimes, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('State Size (nx)');
ylabel('Average Runtime (seconds)');
title('Runtime vs State Size');
grid on;
set(gca, 'XScale', 'log', 'YScale', 'log');

% Subplot 4: Success rate comparison
subplot(2,2,4);
cov_success_rates = zeros(length(state_sizes), 1);
qr_success_rates = zeros(length(state_sizes), 1);

for i = 1:length(state_sizes)
    cov_success_rates(i) = sum(strcmp(results.covariance_steering_status(i, :), 'Success')) / num_trials * 100;
    qr_success_rates(i) = sum(strcmp(results.sqrt_qr_status(i, :), 'Solved')) / num_trials * 100;
end

x = 1:length(state_sizes);
width = 0.35;
bar(x - width/2, cov_success_rates, width, 'FaceColor', [0.2, 0.6, 0.8], 'DisplayName', 'CovarianceSteering');
hold on;
bar(x + width/2, qr_success_rates, width, 'FaceColor', [0.8, 0.4, 0.2], 'DisplayName', 'SqrtQRCovarianceSteering');
xlabel('State Size (nx)');
ylabel('Success Rate (%)');
title('Success Rate Comparison');
legend('Location', 'best');
grid on;
set(gca, 'XTick', 1:length(state_sizes), 'XTickLabel', state_sizes);
ylim([0, 105]);

% Create detailed status table as a separate figure
figure('Position', [200, 200, 1200, 400]);
% Create a table for better visualization
table_data = {};
for i = 1:length(state_sizes)
    for trial = 1:num_trials
        % Format optimal values
        cov_opt_str = 'N/A';
        if ~isnan(results.covariance_steering_optimal(i, trial))
            cov_opt_str = sprintf('%.2e', results.covariance_steering_optimal(i, trial));
        end
        
        sqrt_qr_opt_str = 'N/A';
        if ~isnan(results.sqrt_qr_optimal(i, trial))
            sqrt_qr_opt_str = sprintf('%.2e', results.sqrt_qr_optimal(i, trial));
        end
        
        table_data{end+1, 1} = sprintf('nx=%d, Trial %d', state_sizes(i), trial);
        table_data{end, 2} = results.covariance_steering_status{i, trial};
        table_data{end, 3} = results.sqrt_qr_status{i, trial};
        table_data{end, 4} = cov_opt_str;
        table_data{end, 5} = sqrt_qr_opt_str;
        table_data{end, 6} = sprintf('%.2f', results.runtimes(i, trial));
    end
end

% Display table
uitable('Data', table_data, ...
    'ColumnName', {'Problem', 'Full Cov Status', 'SqrtQR Status', 'FullCovOptVal', 'SqrtQROptVal', 'Runtime (s)'}, ...
    'Position', [20, 20, 1160, 360], ...
    'ColumnWidth', {150, 150, 150, 150, 150, 100});

% Save results
save('runtime_comparison_results.mat', 'results');
fprintf('\nResults saved to runtime_comparison_results.mat\n');