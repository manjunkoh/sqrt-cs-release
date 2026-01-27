%% Scalability vs Horizon Size: Single Obstacle Planning Problem
% Fixed total horizon time; dt changes with N
% To produce plots without rerunning the simulation, load
% data/horizon_size_scalability_single_obstacle_results.mat
clear; clc;
addpath ./SCvxStar/src/
addpath(genpath('./utils'))
addpath('./src')
addpath('./obstacle_path_planning')
figure_settings

%%
% Fixed total time and horizon sizes to test
T_total = 30.0;                  % total time [s]
N_list = [20, 40, 80, 160];     % different horizon lengths to test
num_trials = 5;                   % trials per N for averaging

% Selective rerun parameters
% Set to empty [] to rerun all, or specify subsets
% Methods: 'blockcholesky', 'fullcov', 'sqrtqr'
% Examples:
%   rerun_methods = {'sqrtqr'};           % Rerun only SqrtQR method
%   rerun_nodes = [40, 80];               % Rerun only N=40 and N=80
%   rerun_trials = [1, 3];                % Rerun only trials 1 and 3
%   rerun_methods = {'blockcholesky', 'sqrtqr'}; rerun_nodes = [40];  % Rerun block and qr for N=40 only
rerun_methods = ["blockcholesky"];  % e.g., {'blockcholesky', 'sqrtqr'} to rerun only these methods
rerun_nodes = [20];     % e.g., [40, 80] to rerun only these N values
rerun_trials = [];    % e.g., [1, 3] to rerun only these trial numbers, or [] for all trials

% Load existing results if available
results_file = 'data/horizon_size_scalability_single_obstacle_results.mat';
results_loaded = false;
loaded_N_list = [];
if exist(results_file, 'file')
    fprintf('Loading existing results from %s...\n', results_file);
    loaded_data = load(results_file);
    if isfield(loaded_data, 'results') && isfield(loaded_data, 'N_list')
        loaded_N_list = loaded_data.N_list;
        % Check if all current N values exist in loaded N_list
        if all(ismember(N_list, loaded_N_list))
            results_loaded = true;
            fprintf('  Existing results found. Loaded N_list: %s, Current N_list: %s\n', ...
                mat2str(loaded_N_list), mat2str(N_list));
            fprintf('  Will use loaded results for overlapping N values.\n');
        else
            fprintf('  Warning: Some N values in current list not found in loaded results.\n');
            fprintf('  Loaded: %s, Current: %s\n', mat2str(loaded_N_list), mat2str(N_list));
            fprintf('  Starting fresh for missing N values...\n');
        end
    end
end

% Problem dimensions for 2D double integrator
nx = 4; % [px, py, vx, vy]
nu = 2; % [ax, ay]
nw = nx;

% Basic parameters
wall_y_pos = 1.5;
obstacle_centers = [5, 0];
obstacle_radii = 1.2;
num_obstacles = 1;
pos_idx = 1:2;

mu_0 = [0; 0; 0; 0];
mu_f = [10; 0; 0; 0];

control_risk = 0.005;
u_max = 0.15;
state_risk = 0.005;

% SCP parameters for the square root method
scp_params = SCPParams();
scp_params.tol_opt = 1E-2;
scp_params.tol_feas = 1E-4;
scp_params.k_max = 100;

% Iterative full covariance parameters
max_iters = 100;
penalty_scalar_x = 100;
penalty_scalar_u = 100;
penalty_increase_ratio = 2;
max_penalty = 1e8;
feasibility_tolerance = 1E-4;
convergence_tolerance = 1E-3;

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
results.fullcov_success = false(length(N_list), num_trials);
results.blockcholesky_success = false(length(N_list), num_trials);
results.deterministic_time = zeros(length(N_list), 1);
results.deterministic_iters = zeros(length(N_list), 1);
results.fullcov_iters = zeros(length(N_list), num_trials);
results.sqrtqr_iters = zeros(length(N_list), num_trials);
% Store covariance trajectories (first successful trial for each N x method)
results.fullcov_P = cell(length(N_list), 1);
results.sqrtqr_P = cell(length(N_list), 1);
results.blockcholesky_P = cell(length(N_list), 1);
results.fullcov_mu = cell(length(N_list), 1);
results.sqrtqr_mu = cell(length(N_list), 1);
results.blockcholesky_mu = cell(length(N_list), 1);
results.fullcov_v = cell(length(N_list), 1);
results.sqrtqr_v = cell(length(N_list), 1);
results.blockcholesky_v = cell(length(N_list), 1);
results.fullcov_P_u = cell(length(N_list), 1);
results.sqrtqr_P_u = cell(length(N_list), 1);
results.blockcholesky_P_u = cell(length(N_list), 1);

% Merge with loaded results if available
if results_loaded
    fprintf('  Merging with existing results...\n');
    loaded_results = loaded_data.results;
    % Create mapping from current N_list indices to loaded N_list indices
    N_map = zeros(length(N_list), 1);
    for i = 1:length(N_list)
        idx_in_loaded = find(loaded_N_list == N_list(i), 1);
        if ~isempty(idx_in_loaded)
            N_map(i) = idx_in_loaded;
        end
    end
    
    % Copy over existing results, mapping indices correctly
    field_names = fieldnames(results);
    loaded_num_trials = size(loaded_results.sqrtqr_time, 2);
    num_trials_to_copy = min(num_trials, loaded_num_trials);
    
    for i = 1:length(field_names)
        field_name = field_names{i};
        if isfield(loaded_results, field_name)
            % For cell arrays - map by N index
            if iscell(results.(field_name)) && iscell(loaded_results.(field_name))
                % Check dimensions match
                loaded_dims = size(loaded_results.(field_name));
                results_dims = size(results.(field_name));
                if length(loaded_dims) == length(results_dims) && loaded_dims(2) == results_dims(2)
                    for iN = 1:length(N_list)
                        if N_map(iN) > 0 && N_map(iN) <= loaded_dims(1)
                            if size(results.(field_name), 2) == 1  % Column vector (N x 1)
                                results.(field_name){iN} = loaded_results.(field_name){N_map(iN)};
                            else  % Matrix (N x trials)
                                for trial = 1:num_trials_to_copy
                                    if trial <= size(loaded_results.(field_name), 2)
                                        results.(field_name){iN, trial} = loaded_results.(field_name){N_map(iN), trial};
                                    end
                                end
                            end
                        end
                    end
                end
            % For numeric/logical arrays - map by N index
            elseif (isnumeric(results.(field_name)) || islogical(results.(field_name))) && ...
                   (isnumeric(loaded_results.(field_name)) || islogical(loaded_results.(field_name)))
                % Check dimensions match
                loaded_dims = size(loaded_results.(field_name));
                results_dims = size(results.(field_name));
                if length(loaded_dims) == length(results_dims) && loaded_dims(2) == results_dims(2)
                    for iN = 1:length(N_list)
                        if N_map(iN) > 0 && N_map(iN) <= loaded_dims(1)
                            if size(results.(field_name), 2) == 1  % Column vector (N x 1)
                                results.(field_name)(iN) = loaded_results.(field_name)(N_map(iN));
                            else  % Matrix (N x trials)
                                for trial = 1:num_trials_to_copy
                                    if trial <= size(loaded_results.(field_name), 2)
                                        results.(field_name)(iN, trial) = loaded_results.(field_name)(N_map(iN), trial);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    fprintf('  Results merged (copied %d trials per N).\n', num_trials_to_copy);
end

% Determine which cases to rerun
% Default: rerun all if rerun_methods/rerun_nodes/rerun_trials are empty
if isempty(rerun_methods)
    rerun_methods = {'blockcholesky', 'fullcov', 'sqrtqr'};
end
if isempty(rerun_nodes)
    rerun_nodes = N_list;
end
if isempty(rerun_trials)
    rerun_trials = 1:num_trials;
end

% Convert method names to lowercase for comparison
rerun_methods = lower(rerun_methods);
should_rerun_block = ismember('blockcholesky', rerun_methods);
should_rerun_full = ismember('fullcov', rerun_methods);
should_rerun_qr = ismember('sqrtqr', rerun_methods);

fprintf('\nScalability across horizon sizes with fixed total time T=%.2f s\n', T_total);
fprintf('Horizon sizes: %s\n', mat2str(N_list));
fprintf('Trials per setting: %d\n', num_trials);
if results_loaded
    fprintf('Rerun settings:\n');
    fprintf('  Methods: %s\n', strjoin(rerun_methods, ', '));
    fprintf('  Nodes: %s\n', mat2str(rerun_nodes));
    fprintf('  Trials: %s\n\n', mat2str(rerun_trials));
else
    fprintf('Running all cases...\n\n');
end

rng(1)

for iN = 1:length(N_list)
    N = N_list(iN);
    dt = T_total / N;
    fprintf('============================================================\n');
    fprintf('Testing N = %d (dt = %.4f s)\n', N, dt);

    % Build discrete-time double integrator for 2D
    I2 = eye(2);
    Z2 = zeros(2);
    A = [I2, dt * I2; Z2, I2];
    B = [0.5 * dt^2 * I2; dt * I2];

    % Time-invariant replication along horizon
    A_sys = repmat(A, [1, 1, N]);
    B_sys = repmat(B, [1, 1, N]);

    % Process noise matrix (scaled with dt for stability of covariance growth)
    q = 0.005;
    G = sqrt(q * dt) * eye(nw);
    G_sys = repmat(G, [1, 1, N]);

    % Covariance boundary conditions
    Sigma_0 = diag([0.1, 0.1, 0.01, 0.01]);
    Sigma_f = diag([0.04, 0.04, 0.01, 0.01]);

    % Solve deterministic problem to get reference trajectory (once per N)
    fprintf('  Solving deterministic problem... ');
    yalmip('clear')
    x = sdpvar(nx, N+1, 'full');
    u = sdpvar(nu, N, 'full');
    slack = sdpvar(1, N, 'full');

    x_ref = linspace_vec(mu_0, mu_f, N+1);
    slack_penalty_obstacle = 100;
    max_iters_det = 50;
    flag_deterministic = false;

    t_det = tic;
    for iter = 1:max_iters_det
        constraints = [];
        objective = 0;

        for k = 1:N
            constraints = [constraints, x(:,k+1) == A * x(:,k) + B * u(:,k)];
            constraints = [constraints, u(:,k) <= u_max];
            constraints = [constraints, u(:,k) >= -u_max];
        end

        constraints = [constraints, x(:,1) == mu_0, x(:,N+1) == mu_f];
        constraints = [constraints, x(2,:) <= wall_y_pos];

        for k = 1:N
            for i = 1:size(obstacle_centers, 1)
                a = - (x_ref(pos_idx,k) - obstacle_centers(i,:)');
                b = - 0.5 * norm(x_ref(pos_idx,k) - obstacle_centers(i,:)')^2  + 0.5 * obstacle_radii(i)^2 - a' * x_ref(pos_idx,k);
                constraints = [constraints
                    a' * x(pos_idx,k) + b <= slack(k)
                ];
            end
        end

        constraints = [constraints, slack >= 0];

        for k = 1:N
            objective = objective + norm(u(:,k)) + slack_penalty_obstacle * slack(k);
        end

        sol = optimize(constraints, objective, sdpsettings('verbose', 0));

        if sol.problem
            break;
        end

        x_opt = value(x);
        u_opt = value(u);

        if norm(x_opt - x_ref) < 1e-3 && is_collision_free(x_opt, obstacle_centers, obstacle_radii)
            flag_deterministic = true;
            break;
        end

        x_ref = x_opt;
    end
    det_time = toc(t_det);
    results.deterministic_time(iN) = det_time; % Store once per N
    results.deterministic_iters(iN) = iter; % Store once per N

    if ~flag_deterministic
        fprintf('Deterministic problem failed, skipping N=%d\n', N);
        continue;
    end
    fprintf('done (%.2fs)\n', det_time);

    % Generate hyperplane constraints for each node&obstacle based on the solution
    chance_constraints_state = {};
    for k = 1:N
        for i = 1:num_obstacles
            [a, b] = hyperplane_from_circular_obstacle(obstacle_centers(i,:)', obstacle_radii(i), x_opt(pos_idx,k));
            a = [a; 0; 0];
            b = -b;
            chance_constraints_state = [chance_constraints_state, struct('type', 'affine', 'alpha', a, 'beta', b, 'p', state_risk, 'nodes', k)];
        end
    end

    % Define the chance constraints for control
    chance_constraints_control = {...
        struct('type', 'affine', 'alpha', [1; 0], 'beta', u_max, 'p', control_risk, 'nodes', 1:N), ...
        struct('type', 'affine', 'alpha', [-1; 0], 'beta', u_max, 'p', control_risk, 'nodes', 1:N), ...
        struct('type', 'affine', 'alpha', [0; 1], 'beta', u_max, 'p', control_risk, 'nodes', 1:N), ...
        struct('type', 'affine', 'alpha', [0; -1], 'beta', u_max, 'p', control_risk, 'nodes', 1:N)};

    for trial = 1:num_trials
        % Check if we should rerun this case
        should_rerun_this_case = ismember(N, rerun_nodes) && ismember(trial, rerun_trials);
        
        % Track which methods were actually run
        ran_block = false;
        ran_full = false;
        ran_qr = false;
        
        fprintf('  Trial %2d/%2d ... ', trial, num_trials);

        % Step 2: Solve with Block Cholesky method (no iterations needed)
        if should_rerun_block && should_rerun_this_case
            ran_block = true;
        prob_block = BlockCholeskySteering(...
            A=A_sys, B=B_sys, G=G_sys, ...
            P_0=Sigma_0, P_f=Sigma_f, ...
            objective_type='DV99', ...
            mu_0=mu_0, mu_f=mu_f, ...
            chance_constraints_state=chance_constraints_state, ...
            chance_constraints_control=chance_constraints_control, ...
            N=N);

        t_block = tic;
        diag_block = prob_block.solve(sdpsettings('verbose', 0));
        results.blockcholesky_time(iN, trial) = toc(t_block);

        % Interpret status and optimal value
        switch diag_block.problem
            case 0
                results.blockcholesky_status{iN, trial} = 'SOLVED';
                results.blockcholesky_opt(iN, trial) = prob_block.optimal_objective / N;
                results.blockcholesky_success(iN, trial) = true;
                % Store results from first successful trial
                if isempty(results.blockcholesky_P{iN})
                    results.blockcholesky_P{iN} = prob_block.P;
                    results.blockcholesky_mu{iN} = prob_block.mu;
                    results.blockcholesky_v{iN} = prob_block.v;
                    results.blockcholesky_P_u{iN} = prob_block.P_u;
                end
            case 1
                results.blockcholesky_status{iN, trial} = 'INFEASIBLE';
            case 4
                results.blockcholesky_status{iN, trial} = 'NUMERICAL';
                results.blockcholesky_opt(iN, trial) = prob_block.optimal_objective / N;
                results.blockcholesky_success(iN, trial) = true;
                % Store results even if numerical issues
                if isempty(results.blockcholesky_P{iN})
                    results.blockcholesky_P{iN} = prob_block.P;
                    results.blockcholesky_mu{iN} = prob_block.mu;
                    results.blockcholesky_v{iN} = prob_block.v;
                    results.blockcholesky_P_u{iN} = prob_block.P_u;
                end
            otherwise
                results.blockcholesky_status{iN, trial} = sprintf('ERR%d', diag_block.problem);
        end
        end

        % Step 3: Solve with iterative full covariance method
        if should_rerun_full && should_rerun_this_case
            ran_full = true;
        yalmip('clear')
        mu = sdpvar(nx, N+1, 'full');
        v = sdpvar(nu, N, 'full');
        P = sdpvar(nx, nx, N+1);
        U = sdpvar(nu, nx, N, 'full');
        Y = sdpvar(nu, nu, N);
        slack_u = sdpvar(4, N, 'full');
        slack_x = sdpvar(num_obstacles, N+1, 'full');
        slack_J = sdpvar(1, N);

        mu_ref = x_opt;
        v_ref = u_opt;
        slack_J_ref = (u_max*0.1)^2 * ones(1, N);

        objective_full_covariance = NaN;
        time_full_covariance = NaN;
        flag_fullcov = false;

        penalty_scalar_x_prev = NaN;
        penalty_scalar_u_prev = NaN;
        objective_augmented_prev = Inf;

        t_full = tic;
        for iter = 1:max_iters
            constraints = [];
            objective = 0;

            for k = 1:N
                constraints = [constraints, mu(:,k+1) == A * mu(:,k) + B * v(:,k)];
                constraints = [constraints, P(:,:,k+1) == A * P(:,:,k) * A' + B * Y(:,:,k) * B' + A * U(:,:,k)' * B' + B * U(:,:,k) * A' + G * G'];
                constraints = [constraints, [P(:,:,k) , U(:,:,k)';
                                            U(:,:,k) , Y(:,:,k)] >= 0];
                constraints = [constraints, Y(:,:,k) >= 0];
            end

            constraints = [constraints, mu(:,1) == mu_0, mu(:,N+1) == mu_f, P(:,:,1) == Sigma_0, P(:,:,N+1) <= Sigma_f];

            % Control chance constraints
            for k = 1:N
                for i = 1:length(chance_constraints_control)
                    a = chance_constraints_control{i}.alpha;
                    b = chance_constraints_control{i}.beta;
                    p = chance_constraints_control{i}.p;
                    z = norminv(1 - p);
                    constraints = [constraints
                        z^2 * (a' * Y(:,:,k) * a) <= (b - a' * v_ref(:,k))^2 - 2 * (b - a'* v_ref(:,k)) * a' * (v(:,k) - v_ref(:,k)) + slack_u(i,k)
                        b - a' * v(:,k) >= 0
                    ];
                end
            end

            z = norminv(1 - state_risk);
            for i = 1:num_obstacles
                for k = 1:N
                    [a, b] = hyperplane_from_circular_obstacle(obstacle_centers(i,:)', obstacle_radii(i), x_opt(pos_idx,k));
                    a = [a; 0; 0];
                    b = -b;

                    constraints = [constraints
                        z^2 * (a' * P(:,:,k) * a) <= (b - a'*mu_ref(:,k))^2 - 2 * (b - a'* mu_ref(:,k)) * a' * (mu(:,k) - mu_ref(:,k)) + slack_x(i,k)
                        b - a' * mu(:,k) >= 0
                    ];
                end
            end

            for k = 1:N
                constraints = [constraints
                    lambda_max(Y(:,:,k)) <= slack_J_ref(k)^2 + 2 * slack_J_ref(k) * (slack_J(k) - slack_J_ref(k))
                ];
            end

            constraints = [constraints, slack_u(:) >= 0, slack_x(:) >= 0, slack_J >= 0];

            for k = 1:N
                objective = objective + norm(v(:,k)) + sqrt(chi2inv(0.99, nu)) * slack_J(k) + 1E3 * trace(Y(:,:,k));
            end

            objective_augmented = objective + penalty_scalar_u * sum(slack_u, 'all') + penalty_scalar_x * sum(slack_x, 'all');

            sol = optimize(constraints, objective_augmented, sdpsettings('verbose', 0));

            if sol.problem
                time_full_covariance = toc(t_full);
                results.fullcov_status{iN, trial} = 'INFEASIBLE';
                results.fullcov_iters(iN, trial) = iter;
                % Don't store results if infeasible
                break;
            end

            mu_opt = value(mu);
            v_opt = value(v);
            P_opt = value(P);
            Y_opt = value(Y);
            slack_u_opt = value(slack_u);
            slack_x_opt = value(slack_x);
            
            objective_full_covariance = value(objective);
            objective_augmented = value(objective_augmented);
            objective_with_previous_penalty = objective_full_covariance ...
                + penalty_scalar_x_prev * sum(slack_x_opt, 'all') ...
                + penalty_scalar_u_prev * sum(slack_u_opt, 'all');

            control_constraint_values = check_control_constraint_satisfaction(v_opt, Y_opt, chance_constraints_control, N);
            state_constraint_values = check_state_constraint_satisfaction(mu_opt, P_opt, obstacle_centers, obstacle_radii, x_opt, pos_idx, N, num_obstacles, state_risk);

            if iter > 1 ...
                && (all(slack_x_opt(:) <= feasibility_tolerance) ...
                && all(slack_u_opt(:) <= feasibility_tolerance)) ...
                && all(control_constraint_values(:) <= feasibility_tolerance) ...
                && all(state_constraint_values(:) <= feasibility_tolerance)...
                && objective_augmented_prev - objective_with_previous_penalty < convergence_tolerance
                
                time_full_covariance = toc(t_full);
                % Recalculate objective value after convergence (using actual lambda_max)
                objective_full_covariance_recalc = 0;
                for k = 1:N
                    objective_full_covariance_recalc = objective_full_covariance_recalc + norm(v_opt(:,k)) + sqrt(chi2inv(0.99, nu)) * sqrt(lambda_max(Y_opt(:,:,k)));
                end
                
                results.fullcov_time(iN, trial) = time_full_covariance;
                results.fullcov_opt(iN, trial) = objective_full_covariance_recalc / N;
                results.fullcov_iters(iN, trial) = iter;
                results.fullcov_status{iN, trial} = 'SOLVED';
                results.fullcov_success(iN, trial) = true;
                % Store covariance trajectory from first successful trial
                if isempty(results.fullcov_P{iN})
                    results.fullcov_P{iN} = P_opt;
                    results.fullcov_mu{iN} = mu_opt;
                    results.fullcov_v{iN} = v_opt;
                    results.fullcov_P_u{iN} = Y_opt;
                end
                flag_fullcov = true;
                break;
            end

            v_ref = v_opt;
            mu_ref = mu_opt;

            penalty_scalar_x_prev = penalty_scalar_x;
            penalty_scalar_u_prev = penalty_scalar_u;
            penalty_scalar_x = min(max_penalty, penalty_scalar_x * penalty_increase_ratio);
            penalty_scalar_u = min(max_penalty, penalty_scalar_u * penalty_increase_ratio);
            objective_augmented_prev = objective_augmented;
        end

        if ~flag_fullcov
            if iter == max_iters
                time_full_covariance = toc(t_full);
                results.fullcov_time(iN, trial) = time_full_covariance;
                results.fullcov_iters(iN, trial) = iter;
                results.fullcov_status{iN, trial} = 'MAX_ITERS';
                % Store results even if max iters reached (for plotting)
                % Variables from last iteration should be in scope
                if isempty(results.fullcov_P{iN}) && iter > 0
                    try
                        results.fullcov_P{iN} = P_opt;
                        results.fullcov_mu{iN} = mu_opt;
                        results.fullcov_v{iN} = v_opt;
                        results.fullcov_P_u{iN} = Y_opt;
                    catch
                        % If variables don't exist, skip storing
                    end
                end
            end
        end
        end

        % Step 4: Solve with SqrtQR method
        if should_rerun_qr && should_rerun_this_case
            ran_qr = true;
        init_guess = struct();
        init_guess.S = interpolate_lower_triangular(chol(Sigma_0, 'lower'), chol(Sigma_f, 'lower'), N+1, 'log-cholesky');
        init_guess.L = zeros(nu, nx, N);
        init_guess.mu = x_opt;
        init_guess.v = u_opt;

        prob = SqrtQRCovarianceSteering(init_guess, ...
            N=N, ...
            A_sys=A_sys, B_sys=B_sys, G_sys=G_sys, ...
            objective_type='DV99', ...
            mu_0=mu_0, mu_f=mu_f, ...
            P_0=Sigma_0, P_f=Sigma_f, ...
            chance_constraints_state=chance_constraints_state, ...
            chance_constraints_control=chance_constraints_control);

        flag = prob.solve(save_bool=false, scp_params=scp_params, verbose=false);
        results.sqrtqr_time(iN, trial) = seconds(prob.scp.report.time);
        results.sqrtqr_success(iN, trial) = flag;
        results.sqrtqr_status{iN, trial} = SCPCode.toString(prob.scp.report.code);
        results.sqrtqr_iters(iN, trial) = prob.scp.report.iters;
        % Store results from first trial (successful or not, if postprocess is possible)
        if isempty(results.sqrtqr_P{iN})
            try
                prob.postprocess();
                if ~isempty(prob.P) && ~isempty(prob.mu) && ~isempty(prob.v) && ~isempty(prob.P_u)
                    results.sqrtqr_P{iN} = prob.P;
                    results.sqrtqr_mu{iN} = prob.mu;
                    results.sqrtqr_v{iN} = prob.v;
                    results.sqrtqr_P_u{iN} = prob.P_u;
                end
            catch
                % If postprocess fails, skip storing
            end
        end
        if flag
            results.sqrtqr_opt(iN, trial) = prob.objective(prob.sol) / N;
        end
        end

        % Print summary for this trial
        if ran_block || ran_full || ran_qr
            fprintf('block=%s (%.2fs), full=%s (%.2fs), qr=%s (%.2fs)\n', ...
                results.blockcholesky_status{iN, trial}, results.blockcholesky_time(iN, trial), ...
                results.fullcov_status{iN, trial}, results.fullcov_time(iN, trial), ...
                results.sqrtqr_status{iN, trial}, results.sqrtqr_time(iN, trial));
        else
            fprintf('SKIP (not in rerun list)\n');
        end
    end

    % Summary per N
    valid_qr = results.sqrtqr_time(iN, results.sqrtqr_time(iN, :) > 0);
    valid_full = results.fullcov_time(iN, results.fullcov_time(iN, :) > 0);
    valid_block = results.blockcholesky_time(iN, results.blockcholesky_time(iN, :) > 0);
    if ~isempty(valid_block)
        fprintf('  BlockCholesky avg time: %.2f s\n', mean(valid_block));
    end
    if ~isempty(valid_full)
        fprintf('  FullCov avg time: %.2f s\n', mean(valid_full));
    end
    if ~isempty(valid_qr)
        fprintf('  QR avg time: %.2f s\n\n', mean(valid_qr));
    end
end

%% Console table summary
fprintf('\n=== TERMINATION AND TIME SUMMARY (Fixed T=%.2fs) ===\n', T_total);
fprintf(' N  |   dt    | BlockChol | FullCov | SqrtQR  | BlockCholOptVal | FullCovOptVal | SqrtQROptVal | BlockTime(s) | FullTime(s) | QRTime(s) | FullIters | QRIters\n');
fprintf('----|---------|-----------|---------|---------|-----------------|---------------|--------------|--------------|-------------|-----------|-----------|----------\n');
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
        full_iters_str = 'N/A';
        if results.fullcov_iters(iN, trial) > 0
            full_iters_str = sprintf('%d', results.fullcov_iters(iN, trial));
        end
        qr_iters_str = 'N/A';
        if results.sqrtqr_iters(iN, trial) > 0
            qr_iters_str = sprintf('%d', results.sqrtqr_iters(iN, trial));
        end
        fprintf('%3d | %7.4f | %-9s | %-7s | %-7s | %15s | %13s | %12s | %12.2f | %10.2f | %9.2f | %9s | %8s\n', ...
            N, dt, ...
            results.blockcholesky_status{iN, trial}, ...
            results.fullcov_status{iN, trial}, ...
            results.sqrtqr_status{iN, trial}, ...
            block_opt_str, cov_opt_str, qr_opt_str, ...
            results.blockcholesky_time(iN, trial), ...
            results.fullcov_time(iN, trial), ...
            results.sqrtqr_time(iN, trial), ...
            full_iters_str, qr_iters_str);
    end
    if iN < length(N_list)
        fprintf('----|---------|-----------|---------|---------|-----------------|---------------|--------------|--------------|-------------|-----------|-----------|----------\n');
    end
end


%% Simple plot: average times vs N
avg_full = zeros(length(N_list), 1);
avg_qr = zeros(length(N_list), 1);
avg_block = zeros(length(N_list), 1);
min_full = zeros(length(N_list), 1);
max_full = zeros(length(N_list), 1);
min_qr = zeros(length(N_list), 1);
max_qr = zeros(length(N_list), 1);
min_block = zeros(length(N_list), 1);
max_block = zeros(length(N_list), 1);
for iN = 1:length(N_list)
    % Block Cholesky: min and max
    valid_block = results.blockcholesky_time(iN, results.blockcholesky_time(iN, :) > 0);
    if ~isempty(valid_block)
        avg_block(iN) = mean(valid_block);
        min_block(iN) = min(valid_block);
        max_block(iN) = max(valid_block);
    else
        avg_block(iN) = NaN;
        min_block(iN) = NaN;
        max_block(iN) = NaN;
    end
    
    % Full covariance: min and max
    valid_full = results.fullcov_time(iN, results.fullcov_time(iN, :) > 0);
    if ~isempty(valid_full)
        avg_full(iN) = mean(valid_full);
        min_full(iN) = min(valid_full);
        max_full(iN) = max(valid_full);
    else
        avg_full(iN) = NaN;
        min_full(iN) = NaN;
        max_full(iN) = NaN;
    end
    
    % SqrtQR: average, min and max over valid trials
    times = results.sqrtqr_time(iN, results.sqrtqr_time(iN, :) > 0);
    if ~isempty(times)
        avg_qr(iN) = mean(times);
        min_qr(iN) = min(times);
        max_qr(iN) = max(times);
    else
        avg_qr(iN) = NaN;
        min_qr(iN) = NaN;
        max_qr(iN) = NaN;
    end
end

% Compute average objective values
avg_full_opt = zeros(length(N_list), 1);
avg_qr_opt = zeros(length(N_list), 1);
avg_block_opt = zeros(length(N_list), 1);
for iN = 1:length(N_list)
    % Block Cholesky: average over non-NaN values
    valid_block = results.blockcholesky_opt(iN, ~isnan(results.blockcholesky_opt(iN, :)));
    if ~isempty(valid_block)
        avg_block_opt(iN) = mean(valid_block);
    else
        avg_block_opt(iN) = NaN;
    end
    
    % Full covariance: average over non-NaN values
    valid_full = results.fullcov_opt(iN, ~isnan(results.fullcov_opt(iN, :)));
    if ~isempty(valid_full)
        avg_full_opt(iN) = mean(valid_full);
    else
        avg_full_opt(iN) = NaN;
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
figure(Position=[0, 0, 11, 12]);
hold on;

% Create shaded regions (light grey) for min-max ranges
% Block Cholesky method
valid_block_idx = ~isnan(avg_block);
if any(valid_block_idx)
    x_fill_block = [N_list(valid_block_idx), fliplr(N_list(valid_block_idx))];
    y_fill_block = [min_block(valid_block_idx)', fliplr(max_block(valid_block_idx)')];
    fill(x_fill_block, y_fill_block, "", 'FaceColor', '#D55E00', 'EdgeColor', 'none', 'FaceAlpha', 0.3, HandleVisibility='off');
end

% Full covariance method
valid_full_idx = ~isnan(avg_full);
if any(valid_full_idx)
    x_fill_full = [N_list(valid_full_idx), fliplr(N_list(valid_full_idx))];
    y_fill_full = [min_full(valid_full_idx)', fliplr(max_full(valid_full_idx)')];
    fill(x_fill_full, y_fill_full, "", 'FaceColor', '#0082B2', 'EdgeColor', 'none', 'FaceAlpha', 0.3, HandleVisibility='off');
end

% SqrtQR method (only where valid)
valid_qr_idx = ~isnan(avg_qr);
if any(valid_qr_idx)
    x_fill_qr = [N_list(valid_qr_idx), fliplr(N_list(valid_qr_idx))];
    y_fill_qr = [min_qr(valid_qr_idx)', fliplr(max_qr(valid_qr_idx)')];
    fill(x_fill_qr, y_fill_qr, "", 'FaceColor', '#000000', 'EdgeColor', 'none', 'FaceAlpha', 0.3, HandleVisibility='off');
end

% Plot mean lines
plot(N_list, avg_block, '^-', 'Color', '#D55E00', 'LineWidth', 2, 'MarkerSize', 8, DisplayName='Okamoto \& Tsiotras');
plot(N_list, avg_full, 'o-', 'Color', '#0082B2', 'LineWidth', 2, 'MarkerSize', 8, DisplayName='Liu et al.');
plot(N_list, avg_qr, 's-', 'Color', '#000000', 'LineWidth', 2, 'MarkerSize', 8, DisplayName='Proposed method');

set(gca, 'XScale', 'log', 'YScale', 'log');
grid on;
xlabel('Horizon length N');
ylabel('Average runtime (s)');
xticks(N_list)
legend('Location', 'northwest', 'EdgeColor', 'none', 'BackgroundAlpha', 0.2, 'IconColumnWidth', 15);
exportgraphics(gcf, 'figures/horizon_size_scalability_single_obstacle.png', Resolution=300)
exportgraphics(gcf, 'figures/horizon_size_scalability_single_obstacle.pdf', ContentType='vector')

%% Plot objective function values
figure(Position=[0, 0, 11, 12]);
hold on;

% Plot Block Cholesky method
valid_block_idx = ~isnan(avg_block_opt);
if any(valid_block_idx)
    plot(N_list(valid_block_idx), avg_block_opt(valid_block_idx), '^-', 'Color', '#D55E00', 'LineWidth', 2, 'MarkerSize', 8, DisplayName='Okamoto \& Tsiotras');
end

% Plot full covariance method
valid_full_idx = ~isnan(avg_full_opt);
if any(valid_full_idx)
    plot(N_list(valid_full_idx), avg_full_opt(valid_full_idx), 'o-', 'Color', '#0082B2', 'LineWidth', 2, 'MarkerSize', 8, DisplayName='Liu et al.');
end

% Plot SqrtQR method
valid_qr_idx = ~isnan(avg_qr_opt);
if any(valid_qr_idx)
    plot(N_list(valid_qr_idx), avg_qr_opt(valid_qr_idx), 's-', 'Color', '#000000', 'LineWidth', 2, 'MarkerSize', 8, DisplayName='Proposed method');
end

grid on;
xlabel('Horizon length N');
ylabel('Normalized objective');
xticks(N_list)
legend('Location', 'southwest', 'EdgeColor', 'none', 'BackgroundAlpha', 0.2, 'IconColumnWidth', 15);
xlim([N_list(1) N_list(end)])
xscale log

exportgraphics(gcf, 'figures/horizon_size_cost_comparison_single_obstacle.png', Resolution=300)
exportgraphics(gcf, 'figures/horizon_size_cost_comparison_single_obstacle.pdf', ContentType='vector')

%% Compute average iterations
avg_full_iters = zeros(length(N_list), 1);
avg_qr_iters = zeros(length(N_list), 1);
min_full_iters = zeros(length(N_list), 1);
max_full_iters = zeros(length(N_list), 1);
min_qr_iters = zeros(length(N_list), 1);
max_qr_iters = zeros(length(N_list), 1);
for iN = 1:length(N_list)
    % Full covariance: average, min and max over valid trials
    valid_full_iters = results.fullcov_iters(iN, results.fullcov_iters(iN, :) > 0);
    if ~isempty(valid_full_iters)
        avg_full_iters(iN) = mean(valid_full_iters);
        min_full_iters(iN) = min(valid_full_iters);
        max_full_iters(iN) = max(valid_full_iters);
    else
        avg_full_iters(iN) = NaN;
        min_full_iters(iN) = NaN;
        max_full_iters(iN) = NaN;
    end
    
    % SqrtQR: average, min and max over valid trials
    valid_qr_iters = results.sqrtqr_iters(iN, results.sqrtqr_iters(iN, :) > 0);
    if ~isempty(valid_qr_iters)
        avg_qr_iters(iN) = mean(valid_qr_iters);
        min_qr_iters(iN) = min(valid_qr_iters);
        max_qr_iters(iN) = max(valid_qr_iters);
    else
        avg_qr_iters(iN) = NaN;
        min_qr_iters(iN) = NaN;
        max_qr_iters(iN) = NaN;
    end
end

%% Plot iterations
figure(Position=[0, 0, 11, 12]);
hold on;

% Create shaded regions (light grey) for min-max ranges
% Full covariance method
valid_full_iters_idx = ~isnan(avg_full_iters);
if any(valid_full_iters_idx)
    x_fill_full_iters = [N_list(valid_full_iters_idx), fliplr(N_list(valid_full_iters_idx))];
    y_fill_full_iters = [min_full_iters(valid_full_iters_idx)', fliplr(max_full_iters(valid_full_iters_idx)')];
    fill(x_fill_full_iters, y_fill_full_iters, "", 'FaceColor', '#0082B2', 'EdgeColor', 'none', 'FaceAlpha', 0.3, HandleVisibility='off');
end

% SqrtQR method (only where valid)
valid_qr_iters_idx = ~isnan(avg_qr_iters);
if any(valid_qr_iters_idx)
    x_fill_qr_iters = [N_list(valid_qr_iters_idx), fliplr(N_list(valid_qr_iters_idx))];
    y_fill_qr_iters = [min_qr_iters(valid_qr_iters_idx)', fliplr(max_qr_iters(valid_qr_iters_idx)')];
    fill(x_fill_qr_iters, y_fill_qr_iters, "", 'FaceColor', '#000000', 'EdgeColor', 'none', 'FaceAlpha', 0.3, HandleVisibility='off');
end

% Plot mean lines
plot(N_list, avg_full_iters, 'o-', 'Color', '#0082B2', 'LineWidth', 2, 'MarkerSize', 8, DisplayName='Liu et al.');
plot(N_list, avg_qr_iters, 's-', 'Color', '#000000', 'LineWidth', 2, 'MarkerSize', 8, DisplayName='Proposed method');

set(gca, 'XScale', 'log');
grid on;
xlabel('Horizon length N');
ylabel('Number of iterations');
xticks(N_list)
legend('Location', 'best', 'EdgeColor', 'none', 'BackgroundAlpha', 0.2, 'IconColumnWidth', 15);
xlim([N_list(1) N_list(end)])

exportgraphics(gcf, 'figures/horizon_size_iterations_single_obstacle.png', Resolution=300)
exportgraphics(gcf, 'figures/horizon_size_iterations_single_obstacle.pdf', ContentType='vector')

%% Plot trajectories for each N x method combination
% Store deterministic solution for plotting
results.deterministic_x_opt = cell(length(N_list), 1);
results.deterministic_u_opt = cell(length(N_list), 1);

for iN = 1:length(N_list)
    N = N_list(iN);
    
    % Recompute deterministic solution for this N (needed for plotting)
    dt = T_total / N;
    I2 = eye(2);
    Z2 = zeros(2);
    A = [I2, dt * I2; Z2, I2];
    B = [0.5 * dt^2 * I2; dt * I2];
    
    yalmip('clear')
    x = sdpvar(nx, N+1, 'full');
    u = sdpvar(nu, N, 'full');
    slack = sdpvar(1, N, 'full');
    
    x_ref = linspace_vec(mu_0, mu_f, N+1);
    slack_penalty_obstacle = 100;
    max_iters_det = 50;
    flag_deterministic = false;
    
    for iter = 1:max_iters_det
        constraints = [];
        objective = 0;
        
        for k = 1:N
            constraints = [constraints, x(:,k+1) == A * x(:,k) + B * u(:,k)];
            constraints = [constraints, u(:,k) <= u_max];
            constraints = [constraints, u(:,k) >= -u_max];
        end
        
        constraints = [constraints, x(:,1) == mu_0, x(:,N+1) == mu_f];
        constraints = [constraints, x(2,:) <= wall_y_pos];
        
        for k = 1:N
            for i = 1:size(obstacle_centers, 1)
                a = - (x_ref(pos_idx,k) - obstacle_centers(i,:)');
                b = - 0.5 * norm(x_ref(pos_idx,k) - obstacle_centers(i,:)')^2  + 0.5 * obstacle_radii(i)^2 - a' * x_ref(pos_idx,k);
                constraints = [constraints
                    a' * x(pos_idx,k) + b <= slack(k)
                ];
            end
        end
        
        constraints = [constraints, slack >= 0];
        
        for k = 1:N
            objective = objective + norm(u(:,k)) + slack_penalty_obstacle * slack(k);
        end
        
        sol = optimize(constraints, objective, sdpsettings('verbose', 0));
        
        if sol.problem
            break;
        end
        
        x_opt = value(x);
        u_opt = value(u);
        
        if norm(x_opt - x_ref) < 1e-3 && is_collision_free(x_opt, obstacle_centers, obstacle_radii)
            flag_deterministic = true;
            break;
        end
        
        x_ref = x_opt;
    end
    
    if ~flag_deterministic
        continue;
    end
    
    results.deterministic_x_opt{iN} = x_opt;
    results.deterministic_u_opt{iN} = u_opt;
    
    % Get covariance boundary conditions for this N
    Sigma_0 = diag([0.1, 0.1, 0.01, 0.01]);
    Sigma_f = diag([0.04, 0.04, 0.01, 0.01]);
    
    % Plot for Block Cholesky method
    if ~isempty(results.blockcholesky_P{iN})
        P_block = results.blockcholesky_P{iN};
        mu_block = results.blockcholesky_mu{iN};
        
        figure(Position=[0, 0, 20, 10]);
        hold on;
        plot_obstacles(obstacle_centers, obstacle_radii);
        % Plot initial and final ellipses without legend entries
        plot3sigmaEllipse(mu_0(pos_idx), Sigma_0(pos_idx, pos_idx), Color='#D55E00', HandleVisibility='off');
        plot3sigmaEllipse(mu_f(pos_idx), Sigma_f(pos_idx, pos_idx), Color='#D55E00', LineStyle=":", HandleVisibility='off');
        % Plot covariance ellipses as lines
        for k = 1:size(mu_block, 2)
            P_pos_block = P_block(pos_idx, pos_idx, k);
            if k == 1
                plot3sigmaEllipse(mu_block(pos_idx, k), P_pos_block, 'Color', [.6, .6, .6], DisplayName='$3\sigma$ ellipse');
            else
                plot3sigmaEllipse(mu_block(pos_idx, k), P_pos_block, 'Color', [.6, .6, .6], 'HandleVisibility', 'off');
            end
        end
        % Plot mean trajectory
        plot(mu_block(1,:), mu_block(2,:), 'k.-', DisplayName='Mean', LineWidth=1);
        plot(x_opt(1,:), x_opt(2,:), 'r.-', DisplayName='Reference', LineWidth=1.5);
        % Add text labels at mean locations
        text(mu_0(1), mu_0(2), 'Start', 'FontSize', 20, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
        text(mu_f(1), mu_f(2), 'Goal', 'FontSize', 20, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
        xlabel('$x$');
        ylabel('$y$');
        axis equal;
        set_x_and_y_lims();
        legend(legendUnq(), 'Location', 'northwest', 'EdgeColor', 'none', 'BackgroundAlpha', 0.2, 'IconColumnWidth', 20);
        exportgraphics(gcf, sprintf('figures/trajectory_blockcholesky_N%d.pdf', N), ContentType='vector');
        title(sprintf('Block Cholesky Method - N=%d', N));
        exportgraphics(gcf, sprintf('figures/trajectory_blockcholesky_N%d.png', N), Resolution=300);
    end
    
    % Plot for full covariance method
    if ~isempty(results.fullcov_P{iN})
        P_full = results.fullcov_P{iN};
        mu_full = results.fullcov_mu{iN};
        
        figure(Position=[0, 0, 20, 10]);
        hold on;
        plot_obstacles(obstacle_centers, obstacle_radii);
        % Plot initial and final ellipses without legend entries
        plot3sigmaEllipse(mu_0(pos_idx), Sigma_0(pos_idx, pos_idx), Color='#D55E00', HandleVisibility='off');
        plot3sigmaEllipse(mu_f(pos_idx), Sigma_f(pos_idx, pos_idx), Color='#D55E00', LineStyle=":", HandleVisibility='off');
        % Plot covariance ellipses as lines
        for k = 1:size(mu_full, 2)
            P_pos_full = P_full(pos_idx, pos_idx, k);
            if k == 1
                plot3sigmaEllipse(mu_full(pos_idx, k), P_pos_full, 'Color', [.6, .6, .6], DisplayName='$3\sigma$ ellipse');
            else
                plot3sigmaEllipse(mu_full(pos_idx, k), P_pos_full, 'Color', [.6, .6, .6], 'HandleVisibility', 'off');
            end
        end
        % Plot mean trajectory
        plot(mu_full(1,:), mu_full(2,:), 'k.-', DisplayName='Mean', LineWidth=1);
        plot(x_opt(1,:), x_opt(2,:), 'r.-', DisplayName='Reference', LineWidth=1.5);
        % Add text labels at mean locations
        text(mu_0(1), mu_0(2), 'Start', 'FontSize', 20, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
        text(mu_f(1), mu_f(2), 'Goal', 'FontSize', 20, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
        xlabel('$x$');
        ylabel('$y$');
        axis equal;
        set_x_and_y_lims();
        legend(legendUnq(), 'Location', 'northwest', 'EdgeColor', 'none', 'BackgroundAlpha', 0.2, 'IconColumnWidth', 20);
        exportgraphics(gcf, sprintf('figures/trajectory_fullcov_N%d.pdf', N), ContentType='vector');
        title(sprintf('Full Covariance Method - N=%d', N));
        exportgraphics(gcf, sprintf('figures/trajectory_fullcov_N%d.png', N), Resolution=300);
    end
    
    % Plot for SqrtQR method
    if ~isempty(results.sqrtqr_P{iN})
        P_qr = results.sqrtqr_P{iN};
        mu_qr = results.sqrtqr_mu{iN};
        
        figure(Position=[0, 0, 20, 10]);
        hold on;
        plot_obstacles(obstacle_centers, obstacle_radii);
        % Plot initial and final ellipses without legend entries
        plot3sigmaEllipse(mu_0(pos_idx), Sigma_0(pos_idx, pos_idx), Color='#D55E00', HandleVisibility='off');
        plot3sigmaEllipse(mu_f(pos_idx), Sigma_f(pos_idx, pos_idx), Color='#D55E00', LineStyle=":", HandleVisibility='off');
        % Plot covariance ellipses as lines
        for k = 1:size(mu_qr, 2)
            P_pos_qr = P_qr(pos_idx, pos_idx, k);
            if k == 1
                plot3sigmaEllipse(mu_qr(pos_idx, k), P_pos_qr, 'Color', [.6, .6, .6], DisplayName='$3\sigma$ ellipse');
            else
                plot3sigmaEllipse(mu_qr(pos_idx, k), P_pos_qr, 'Color', [.6, .6, .6], 'HandleVisibility', 'off');
            end
        end
        % Plot mean trajectory
        plot(mu_qr(1,:), mu_qr(2,:), 'k.-', DisplayName='Mean', LineWidth=1);
        plot(x_opt(1,:), x_opt(2,:), 'r.-', DisplayName='Reference', LineWidth=1.5);
        % Add text labels at mean locations
        text(mu_0(1), mu_0(2), 'Start', 'FontSize', 20, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
        text(mu_f(1), mu_f(2), 'Goal', 'FontSize', 20, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
        xlabel('$x$');
        ylabel('$y$');
        axis equal;
        set_x_and_y_lims();
        legend(legendUnq(), 'Location', 'northwest', 'EdgeColor', 'none', 'BackgroundAlpha', 0.2, 'IconColumnWidth', 20);
        exportgraphics(gcf, sprintf('figures/trajectory_sqrtqr_N%d.pdf', N), ContentType='vector');
        title(sprintf('SqrtQR Method - N=%d', N));
        exportgraphics(gcf, sprintf('figures/trajectory_sqrtqr_N%d.png', N), Resolution=300);
    end
end

function set_x_and_y_lims()
    xlim([-1, 11]);
    ylim([-1, 5]);
end

%% Plot control trajectories for each N x method combination
for iN = 1:length(N_list)
    N = N_list(iN);
    dt = T_total / N;
    t_vec = (0:N-1) * dt;  % Control is at nodes 0 to N-1
    
    % Plot for Block Cholesky method
    if ~isempty(results.blockcholesky_v{iN}) && ~isempty(results.blockcholesky_P_u{iN}) && ~isempty(results.deterministic_u_opt{iN})
        v_block = results.blockcholesky_v{iN};
        P_u_block = results.blockcholesky_P_u{iN};
        u_ref = results.deterministic_u_opt{iN};
        
        % Compute standard deviations for each control component
        sigma_u1_block = zeros(1, N);
        sigma_u2_block = zeros(1, N);
        for k = 1:N
            sigma_u1_block(k) = sqrt(P_u_block(1, 1, k));
            sigma_u2_block(k) = sqrt(P_u_block(2, 2, k));
        end
        
        % Compute 99% bounds for each component using norminv
        z_99 = norminv(0.99);
        bound99_u1_upper = v_block(1,:) + z_99 * sigma_u1_block;
        bound99_u1_lower = v_block(1,:) - z_99 * sigma_u1_block;
        bound99_u2_upper = v_block(2,:) + z_99 * sigma_u2_block;
        bound99_u2_lower = v_block(2,:) - z_99 * sigma_u2_block;
        
        % Compute 99% bound for norm
        q_99 = sqrt(chi2inv(0.99, nu));
        norm_bound_block = zeros(1, N);
        norm_mean_block = zeros(1, N);
        for k = 1:N
            norm_mean_block(k) = norm(v_block(:, k));
            norm_bound_block(k) = norm_mean_block(k) + q_99 * sqrt(lambda_max(P_u_block(:,:,k)));
        end
        
        % Compute reference norm
        norm_ref = zeros(1, N);
        for k = 1:N
            norm_ref(k) = norm(u_ref(:, k));
        end
        
        figure(Position=[0, 0, 18, 18]);
        tiledlayout(3,1);
        
        % Plot u1 (x-direction)
        nexttile;
        hold on;
        plot(t_vec, v_block(1,:), '^-', 'Color', '#D55E00', 'LineWidth', 2, 'MarkerSize', 6, DisplayName='Mean');
        plot(t_vec, bound99_u1_upper, '--', 'Color', '#D55E00', 'LineWidth', 1.5, DisplayName='99% bound');
        plot(t_vec, bound99_u1_lower, '--', 'Color', '#D55E00', 'LineWidth', 1.5, HandleVisibility='off');
        plot(t_vec, u_ref(1,:), 'r-', 'LineWidth', 1.5, DisplayName='Reference');
        yline(u_max, 'k:', 'LineWidth', 1.5, DisplayName='$u_{\max}$');
        yline(-u_max, 'k:', 'LineWidth', 1.5, DisplayName='$-u_{\max}$');
        grid on;
        ylabel('Control $u_x$');
        legend(legendUnq(), 'Location', 'northeast', 'Orientation', 'horizontal', 'EdgeColor', 'none', 'IconColumnWidth', 15);
        
        % Plot u2 (y-direction)
        nexttile;
        hold on;
        plot(t_vec, v_block(2,:), '^-', 'Color', '#D55E00', 'LineWidth', 2, 'MarkerSize', 6, DisplayName='Mean');
        plot(t_vec, bound99_u2_upper, '--', 'Color', '#D55E00', 'LineWidth', 1.5, DisplayName='99% bound');
        plot(t_vec, bound99_u2_lower, '--', 'Color', '#D55E00', 'LineWidth', 1.5, HandleVisibility='off');
        plot(t_vec, u_ref(2,:), 'r-', 'LineWidth', 1.5, DisplayName='Reference');
        yline(u_max, 'k:', 'LineWidth', 1.5, DisplayName='$u_{\max}$');
        yline(-u_max, 'k:', 'LineWidth', 1.5, DisplayName='$-u_{\max}$');
        grid on;
        ylabel('Control $u_y$');

        % Plot norm with 99% bound
        nexttile;
        hold on;
        plot(t_vec, norm_mean_block, '^-', 'Color', '#D55E00', 'LineWidth', 2, 'MarkerSize', 6, DisplayName='Mean');
        plot(t_vec, norm_bound_block, '--', 'Color', '#D55E00', 'LineWidth', 1.5, DisplayName='99% bound');
        plot(t_vec, norm_ref, 'r-', 'LineWidth', 1.5, DisplayName='Reference');
        grid on;
        xlabel('Time (s)');
        ylabel('$\|u\|_2$');
        ylim([0, 0.25]);

        exportgraphics(gcf, sprintf('figures/control_blockcholesky_N%d.pdf', N), ContentType='vector');        
        sgtitle(sprintf('Control Trajectory - Block Cholesky Method, N=%d', N));
        exportgraphics(gcf, sprintf('figures/control_blockcholesky_N%d.png', N), Resolution=300);
    end
    
    % Plot for full covariance method
    if ~isempty(results.fullcov_v{iN}) && ~isempty(results.fullcov_P_u{iN}) && ~isempty(results.deterministic_u_opt{iN})
        v_full = results.fullcov_v{iN};
        P_u_full = results.fullcov_P_u{iN};
        u_ref = results.deterministic_u_opt{iN};
        
        % Compute standard deviations for each control component
        sigma_u1_full = zeros(1, N);
        sigma_u2_full = zeros(1, N);
        for k = 1:N
            sigma_u1_full(k) = sqrt(P_u_full(1, 1, k));
            sigma_u2_full(k) = sqrt(P_u_full(2, 2, k));
        end
        
        % Compute 99% bounds for each component using norminv
        z_99 = norminv(0.99);
        bound99_u1_upper = v_full(1,:) + z_99 * sigma_u1_full;
        bound99_u1_lower = v_full(1,:) - z_99 * sigma_u1_full;
        bound99_u2_upper = v_full(2,:) + z_99 * sigma_u2_full;
        bound99_u2_lower = v_full(2,:) - z_99 * sigma_u2_full;
        
        % Compute 99% bound for norm
        q_99 = sqrt(chi2inv(0.99, nu));
        norm_bound_full = zeros(1, N);
        norm_mean_full = zeros(1, N);
        for k = 1:N
            norm_mean_full(k) = norm(v_full(:, k));
            norm_bound_full(k) = norm_mean_full(k) + q_99 * sqrt(lambda_max(P_u_full(:,:,k)));
        end
        
        % Compute reference norm
        norm_ref = zeros(1, N);
        for k = 1:N
            norm_ref(k) = norm(u_ref(:, k));
        end
        
        figure(Position=[0, 0, 18, 18]);
        tiledlayout(3,1);
        
        % Plot u1 (x-direction)
        nexttile;
        hold on;
        plot(t_vec, v_full(1,:), 'o-', 'Color', '#0082B2', 'LineWidth', 2, 'MarkerSize', 6, DisplayName='Mean');
        plot(t_vec, bound99_u1_upper, '--', 'Color', '#0082B2', 'LineWidth', 1.5, DisplayName='99\% bound');
        plot(t_vec, bound99_u1_lower, '--', 'Color', '#0082B2', 'LineWidth', 1.5, HandleVisibility='off');
        plot(t_vec, u_ref(1,:), 'r-', 'LineWidth', 1.5, DisplayName='Reference');
        yline(u_max, 'k:', 'LineWidth', 1.5, HandleVisibility='off');
        yline(-u_max, 'k:', 'LineWidth', 1.5, HandleVisibility='off');
        grid on;
        ylabel('Control $u_x$');
        legend(legendUnq(), 'Location', 'northeast', 'Orientation', 'horizontal', 'EdgeColor', 'none', 'IconColumnWidth', 15);
        
        % Plot u2 (y-direction)
        nexttile;
        hold on;
        plot(t_vec, v_full(2,:), 'o-', 'Color', '#0082B2', 'LineWidth', 2, 'MarkerSize', 6, DisplayName='Mean');
        plot(t_vec, bound99_u2_upper, '--', 'Color', '#0082B2', 'LineWidth', 1.5, DisplayName='99\% bound');
        plot(t_vec, bound99_u2_lower, '--', 'Color', '#0082B2', 'LineWidth', 1.5, HandleVisibility='off');
        plot(t_vec, u_ref(2,:), 'r-', 'LineWidth', 1.5, DisplayName='Reference');
        yline(u_max, 'k:', 'LineWidth', 1.5, HandleVisibility='off');
        yline(-u_max, 'k:', 'LineWidth', 1.5, HandleVisibility='off');
        grid on;
        ylabel('Control $u_y$');
        
        % Plot norm with 99% bound
        nexttile;
        hold on;
        plot(t_vec, norm_mean_full, 'o-', 'Color', '#0082B2', 'LineWidth', 2, 'MarkerSize', 6, DisplayName='Mean');
        plot(t_vec, norm_bound_full, '--', 'Color', '#0082B2', 'LineWidth', 1.5, DisplayName='99\% bound');
        plot(t_vec, norm_ref, 'r-', 'LineWidth', 1.5, DisplayName='Reference');
        grid on;
        xlabel('Time (s)');
        ylabel('$\|u\|_2$');
        ylim([0, 0.25]);

        exportgraphics(gcf, sprintf('figures/control_fullcov_N%d.pdf', N), ContentType='vector');        
        sgtitle(sprintf('Control Trajectory - Full Covariance Method, N=%d', N));
        exportgraphics(gcf, sprintf('figures/control_fullcov_N%d.png', N), Resolution=300);
    end
    
    % Plot for SqrtQR method
    if ~isempty(results.sqrtqr_v{iN}) && ~isempty(results.sqrtqr_P_u{iN}) && ~isempty(results.deterministic_u_opt{iN})
        v_qr = results.sqrtqr_v{iN};
        P_u_qr = results.sqrtqr_P_u{iN};
        u_ref = results.deterministic_u_opt{iN};
        
        % Compute standard deviations for each control component
        sigma_u1_qr = zeros(1, N);
        sigma_u2_qr = zeros(1, N);
        for k = 1:N
            sigma_u1_qr(k) = sqrt(P_u_qr(1, 1, k));
            sigma_u2_qr(k) = sqrt(P_u_qr(2, 2, k));
        end
        
        % Compute 99% bounds for each component using norminv
        z_99 = norminv(0.99);
        bound99_u1_upper = v_qr(1,:) + z_99 * sigma_u1_qr;
        bound99_u1_lower = v_qr(1,:) - z_99 * sigma_u1_qr;
        bound99_u2_upper = v_qr(2,:) + z_99 * sigma_u2_qr;
        bound99_u2_lower = v_qr(2,:) - z_99 * sigma_u2_qr;
        
        % Compute 99% bound for norm
        q_99 = sqrt(chi2inv(0.99, nu));
        norm_bound_qr = zeros(1, N);
        norm_mean_qr = zeros(1, N);
        for k = 1:N
            norm_mean_qr(k) = norm(v_qr(:, k));
            norm_bound_qr(k) = norm_mean_qr(k) + q_99 * sqrt(lambda_max(P_u_qr(:,:,k)));
        end
        
        % Compute reference norm
        norm_ref = zeros(1, N);
        for k = 1:N
            norm_ref(k) = norm(u_ref(:, k));
        end
        
        figure(Position=[0, 0, 18, 18]);
        tiledlayout(3,1);
        
        % Plot u1 (x-direction)
        nexttile;
        hold on;
        plot(t_vec, v_qr(1,:), 's-', 'Color', '#000000', 'LineWidth', 2, 'MarkerSize', 6, DisplayName='Mean');
        plot(t_vec, bound99_u1_upper, '--', 'Color', '#000000', 'LineWidth', 1.5, DisplayName='99\% bound');
        plot(t_vec, bound99_u1_lower, '--', 'Color', '#000000', 'LineWidth', 1.5, HandleVisibility='off');
        plot(t_vec, u_ref(1,:), 'r-', 'LineWidth', 1.5, DisplayName='Reference');
        yline(u_max, 'k:', 'LineWidth', 1.5, HandleVisibility='off');
        yline(-u_max, 'k:', 'LineWidth', 1.5, HandleVisibility='off');
        grid on;
        ylabel('Control $u_x$');
        legend(legendUnq(), 'Location', 'northeast', 'Orientation', 'horizontal', 'EdgeColor', 'none', 'IconColumnWidth', 15);
        
        % Plot u2 (y-direction)
        nexttile;
        hold on;
        plot(t_vec, v_qr(2,:), 's-', 'Color', '#000000', 'LineWidth', 2, 'MarkerSize', 6, DisplayName='Mean');
        plot(t_vec, bound99_u2_upper, '--', 'Color', '#000000', 'LineWidth', 1.5, DisplayName='99\% bound');
        plot(t_vec, bound99_u2_lower, '--', 'Color', '#000000', 'LineWidth', 1.5, HandleVisibility='off');
        plot(t_vec, u_ref(2,:), 'r-', 'LineWidth', 1.5, DisplayName='Reference');
        yline(u_max, 'k:', 'LineWidth', 1.5, HandleVisibility='off');
        yline(-u_max, 'k:', 'LineWidth', 1.5, HandleVisibility='off');
        grid on;
        ylabel('Control $u_y$');

        % Plot norm with 99% bound
        nexttile;
        hold on;
        plot(t_vec, norm_mean_qr, 's-', 'Color', '#000000', 'LineWidth', 2, 'MarkerSize', 6, DisplayName='Mean');
        plot(t_vec, norm_bound_qr, '--', 'Color', '#000000', 'LineWidth', 1.5, DisplayName='99% bound');
        plot(t_vec, norm_ref, 'r-', 'LineWidth', 1.5, DisplayName='Reference');
        grid on;
        xlabel('Time (s)');
        ylabel('$\|u\|_2$');
        ylim([0, 0.25]);
        
        exportgraphics(gcf, sprintf('figures/control_sqrtqr_N%d.pdf', N), ContentType='vector');        
        sgtitle(sprintf('Control Trajectory - SqrtQR Method, N=%d', N));
        exportgraphics(gcf, sprintf('figures/control_sqrtqr_N%d.png', N), Resolution=300);
    end
end

%%
% Create data directory if it doesn't exist
if ~exist('data', 'dir')
    mkdir('data')
end
save('data/horizon_size_scalability_single_obstacle_results.mat', 'results', 'T_total', 'N_list');
fprintf('\nSaved results to data/horizon_size_scalability_single_obstacle_results.mat\n');
