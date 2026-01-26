clc; clear;
addpath(genpath('./utils'))
addpath ./src
addpath ./SCvxStar/src/
addpath ./obstacle_path_planning
figure_settings

%% Parameters
basic_parameters;

% Define the chance constraints
chance_constraints_control={...
    struct('type', 'affine', 'alpha', [0; 1], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes), ...
    struct('type', 'affine', 'alpha', [0; -1], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes), ...
    struct('type', 'affine', 'alpha', [1; 0], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes), ...
    struct('type', 'affine', 'alpha', [-1; 0], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes)};

chance_constraints_state={struct('type', 'affine', 'alpha', ...
    [0; 1; 0; 0], 'beta', wall_y_pos, 'p', state_risk, 'nodes', 1:num_nodes+1)};

num_simulations = 10000;
%% Load generated obstacle environments
load_filename = './data/obstacle_cases.mat';
if ~exist(load_filename, 'file')
    error('Obstacle cases file not found: %s\nPlease run generate_and_save_obstacles.m first', load_filename);
end

fprintf('Loading obstacle cases from %s...\n', load_filename);
loaded_data = load(load_filename);
valid_cases = loaded_data.valid_cases;
% num_trials = length(valid_cases.obstacle_centers);
num_trials = 5;
num_obstacles = loaded_data.num_obstacles;
fprintf('Running %d trials...\n\n', num_trials);

results = struct();
results.num_simulations = num_simulations;
results.x_opt_deterministic = cell(num_trials, 1);
results.u_opt_deterministic = cell(num_trials, 1);
results.obstacle_centers = cell(num_trials, 1);
results.obstacle_radii = cell(num_trials, 1);

results.iters_full_covariance = NaN(num_trials, 1);
results.objective_full_covariance = NaN(num_trials, 1);
results.time_full_covariance = NaN(num_trials, 1);
results.prob_full_covariance = cell(num_trials, 1);
results.solved_full_covariance = false(num_trials, 1);

results.iters_qr = NaN(num_trials, 1);
results.objective_qr = NaN(num_trials, 1);
results.time_qr = NaN(num_trials, 1);
results.prob_qr = cell(num_trials, 1);
results.solved_qr = false(num_trials, 1);

% Monte Carlo simulation results
results.mc_num_simulations = [];
results.mc_collision_counts_full_covariance = cell(num_trials, 1);
results.mc_collision_counts_qr = cell(num_trials, 1);
results.mc_collision_prob_full_covariance = NaN(num_trials, num_nodes+1, num_obstacles);
results.mc_collision_prob_qr = NaN(num_trials, num_nodes+1, num_obstacles);

%% SCP parameters
scp_params = SCPParams();
scp_params.k_max = 100;
scp_params.tol_opt = 1E-2;
scp_params.tol_feas = 1E-4;
scp_params.linearization = 'inexact';

for trial = 1:num_trials
    fprintf('Trial %d/%d ... ', trial, num_trials);
    
    % Load pre-computed obstacle environment and deterministic solution
    obstacle_centers = valid_cases.obstacle_centers{trial};
    obstacle_radii = valid_cases.obstacle_radii{trial};
    x_opt = valid_cases.x_opt{trial};
    u_opt = valid_cases.u_opt{trial};
    
    % Ensure obstacle_radii is a row vector for consistent indexing
    if size(obstacle_radii, 1) > size(obstacle_radii, 2)
        obstacle_radii = obstacle_radii';
    end
    
    % Store in results
    results.obstacle_centers{trial} = obstacle_centers;
    results.obstacle_radii{trial} = obstacle_radii;
    results.x_opt_deterministic{trial} = x_opt;
    results.u_opt_deterministic{trial} = u_opt;
        
    % Create circular obstacles structure for stochastic problem
    circular_obstacles = {};
    for i = 1:num_obstacles
        circular_obstacles = [circular_obstacles, struct('center', obstacle_centers(i,:)', 'radius', obstacle_radii(i), 'p', state_risk)];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Solve with Full Covariance formulation
    max_iters = 100;
    penalty_scalar_x = 1000;
    penalty_scalar_u = 100;
    convergence_tolerance = 1E-4;

    yalmip('clear')
    mu = sdpvar(nx, num_nodes+1, 'full');
    v = sdpvar(nu, num_nodes, 'full');
    P = sdpvar(nx, nx, num_nodes+1);
    U = sdpvar(nu, nx, num_nodes, 'full');
    Y = sdpvar(nu, nu, num_nodes);
    lambda = sdpvar(1, num_nodes, 'full'); % "virtual control"
    lambda_x = sdpvar(num_obstacles, num_nodes+1, 'full');

    mu_ref = x_opt;
    v_ref = u_opt;

    Y_init_value = 1e-2; % This value is really important for optimality
    P_ref = interpolate_lower_triangular(chol(Sigma_0, 'lower'), chol(Sigma_f, 'lower'), num_nodes+1, 'log-cholesky');
    Y_ref = repmat(Y_init_value * eye(nu), [1,1,num_nodes]); 

    % Initialize tracking variables
    objective_full_covariance = NaN;
    time_full_covariance = NaN;

    tic;
    for iter = 1:max_iters
        constraints = [];
        objective = 0;

        for k = 1:num_nodes
            constraints = [constraints, mu(:,k+1) == A * mu(:,k) + B * v(:,k)];
            constraints = [constraints, P(:,:,k+1) == A * P(:,:,k) * A' + B * Y(:,:,k) * B' + A * U(:,:,k)' * B' + B * U(:,:,k) * A' + G * G'];
            constraints = [constraints, [P(:,:,k) , U(:,:,k)' ;
                                        U(:,:,k) , Y(:,:,k)] >= 0];
            constraints = [constraints, P(:,:,k) >= 0, Y(:,:,k) >= 0];
        end

        for k = 1:num_nodes + 1
            constraints = [constraints, P(:,:,k) >= 0];
        end

        constraints = [constraints, mu(:,1) == mu_0, mu(:,num_nodes+1) == mu_f, P(:,:,1) == Sigma_0, P(:,:,num_nodes+1) <= Sigma_f];

        z = norminv(1 - state_risk);

        % State wall chance constraints
        for k = 1:num_nodes
            sqrt_ref = sqrt(P_ref(2,2,k));
            constraints = [constraints,
                z / (2 * sqrt_ref) * (P(2,2,k)) + z * sqrt_ref / 2 + mu(2,k) - wall_y_pos <= 0
            ];
        end

        % Control chance constraints
        for k = 1:num_nodes
            v_k = v(:,k);
            Y_k = Y(:,:,k);
            v_ref_k = v_ref(:,k);
            Y_ref_k = Y_ref(:,:,k);
            for i = 1:length(chance_constraints_control)
                a = chance_constraints_control{i}.alpha;
                b = chance_constraints_control{i}.beta;
                p = chance_constraints_control{i}.p;
                z = norminv(1 - p);
                sqrt_ref = sqrt(a' * Y_ref_k * a);
                constraints = [constraints
                    z / (2 * sqrt_ref) * (a' * Y_k * a) + a' * v_k - b + z * sqrt_ref / 2 <= lambda(k)
                ];
            end
        end

        % State obstacle chance constraints
        for k = 1:num_nodes+1
            for i = 1:num_obstacles
                [a, b] = hyperplane_from_circular_obstacle(obstacle_centers(i,:)', obstacle_radii(i), mu_ref(pos_idx,k));
                P_ref_pos_k = P_ref(pos_idx,pos_idx,k);
                sqrt_ref = sqrt(a' * P_ref_pos_k * a);
                if sqrt_ref <= 0
                    % numerically, the covariance can be non-PSD; in this case
                    % simply set this value to a small positive number
                    % Using nearestSPD sometimes does not terminate for a long
                    % time, so it isn't used here
                    sqrt_ref = 0.0001;
                end
                constraints = [constraints
                    z / (2 * sqrt_ref) * (a' * P(pos_idx,pos_idx,k) * a) + z * sqrt_ref / 2 + a' * mu(pos_idx,k) + b <= lambda_x(i,k)
                ];
            end
        end

        constraints = [constraints, lambda >= 0, lambda_x(:) >= 0];

        for k = 1:num_nodes
            objective = objective + mu(:,k)' * Q * mu(:,k) + v(:,k)' * R * v(:,k) + trace(Q * P(:,:,k)) + trace(R * Y(:,:,k));
        end

        objective_augmented = objective + penalty_scalar_u * sum(lambda) + penalty_scalar_x + sum(lambda_x, 'all');

        fprintf("Iteration %d    ", iter)

        sol = optimize(constraints, objective_augmented, sdpsettings('verbose', 0));


        if sol.problem
            time_full_covariance = toc;
            fprintf('Infeasible at iteration %d\n', iter);
            break;
        end

        fprintf("Objective: %f\n",  value(objective))

        mu_opt = value(mu);
        v_opt = value(v);
        P_opt = value(P);
        U_opt = value(U);
        Y_opt = value(Y);
        lambda_opt = value(lambda);
        lambda_x_opt = value(lambda_x);
        
        objective_full_covariance = value(objective);

        if all(vecnorm(mu_opt - mu_ref, Inf) < convergence_tolerance) 
            time_full_covariance = toc;
            fprintf('Converged in %d iterations in %.3f seconds\n', iter, time_full_covariance);
            prob_fc = struct();
            prob_fc.mu = mu_opt;
            prob_fc.v = v_opt;
            prob_fc.P = P_opt;
            prob_fc.Y = Y_opt;
            prob_fc.lambda = lambda_opt;
            prob_fc.lambda_x = lambda_x_opt;
            % Compute feedback gain K for Monte Carlo simulation
            prob_fc.K = zeros(nu, nx, num_nodes);
            for k = 1:num_nodes
                prob_fc.K(:,:,k) = U_opt(:,:,k) / P_opt(:,:,k);
            end
            results.solved_full_covariance(trial) = true;
            break;
        end

        mu_ref = mu_opt;
        P_ref = P_opt;

        if iter == max_iters
            time_full_covariance = toc;
            fprintf('Reached maximum iterations\n')
            % Store solution even if max iterations reached
            if ~exist('prob_fc', 'var') || isempty(prob_fc)
                prob_fc = struct();
                prob_fc.mu = mu_opt;
                prob_fc.v = v_opt;
                prob_fc.P = P_opt;
                prob_fc.Y = Y_opt;
                prob_fc.lambda = lambda_opt;
                prob_fc.lambda_x = lambda_x_opt;
                % Compute feedback gain K for Monte Carlo simulation
                prob_fc.K = zeros(nu, nx, num_nodes);
                for k = 1:num_nodes
                    prob_fc.K(:,:,k) = U_opt(:,:,k) / P_opt(:,:,k);
                end
            end
        end
    end

    iters_full_covariance = iter;
    
    % Initialize prob_fc as empty if it wasn't created
    if ~exist('prob_fc', 'var')
        prob_fc = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve with SQRT QR method

    init_guess = struct();
    init_guess.S = interpolate_lower_triangular(chol(Sigma_0, 'lower'), chol(Sigma_f, 'lower'), num_nodes+1, 'log-cholesky');
    init_guess.L = zeros(nu, nx, num_nodes);
    init_guess.mu = x_opt;
    init_guess.v = u_opt;

    prob_qr = SqrtQRCovarianceSteering(init_guess, ...
        N=num_nodes, ...
        A_sys=A_sys, B_sys=B_sys, G_sys=G_sys, ...
        objective_type='LQR', ...
        mu_0=mu_0, mu_f=mu_f, ...
        P_0=Sigma_0, P_f=Sigma_f, Q=Q, R=R, ...
        chance_constraints_state=chance_constraints_state, ...
        chance_constraints_control=chance_constraints_control,...
        circular_obstacles=circular_obstacles, ...
        relax_obstacle_constraints=false);

    flag_solved_qr = prob_qr.solve(scp_params=scp_params);

    if flag_solved_qr
        prob_qr.postprocess();
        time_qr = seconds(prob_qr.scp.report.time);
        iters_qr = prob_qr.scp.report.iters;
        objective_qr = prob_qr.objective(prob_qr.sol);
        fprintf('SqrtQRCovarianceSteering solved successfully in %.3f seconds\n', time_qr);
        fprintf('  Number of iterations: %d\n', iters_qr);
        fprintf('  Objective value: %.6f\n', objective_qr);
        results.solved_qr(trial) = true;
    else
        time_qr = NaN;
        iters_qr = NaN;
        objective_qr = NaN;
        prob_qr = [];
        fprintf('SqrtQRCovarianceSteering did not solve successfully\n');
    end

    results.iters_full_covariance(trial) = iters_full_covariance;
    results.objective_full_covariance(trial) = objective_full_covariance;
    results.time_full_covariance(trial) = time_full_covariance;
    results.iters_qr(trial) = iters_qr;
    results.objective_qr(trial) = objective_qr;
    results.time_qr(trial) = time_qr;
    results.prob_full_covariance{trial} = prob_fc;
    results.prob_qr{trial} = prob_qr;
    
end

%% Monte Carlo simulation for each method
disp("Running Monte Carlo simulations...")
for trial = 1:num_trials

    prob_fc = results.prob_full_covariance{trial};
    prob_qr = results.prob_qr{trial};

    solved_full_covariance = results.solved_full_covariance(trial);
    solved_qr = results.solved_qr(trial);

    obstacle_centers = results.obstacle_centers{trial};
    obstacle_radii = results.obstacle_radii{trial};

    % Monte Carlo simulation for Full Covariance method
    if solved_full_covariance

        [x_hist_all_fc, u_hist_all_fc] = simulate_samples(mu_0, Sigma_0, A_sys, B_sys, G_sys, prob_fc.K, prob_fc.mu, prob_fc.v, num_nodes, num_simulations);
        collision_counts_fc = count_collision_samples(x_hist_all_fc, obstacle_centers, obstacle_radii);
        results.mc_collision_counts_full_covariance{trial} = collision_counts_fc;
        results.mc_collision_prob_full_covariance(trial, :, :) = collision_counts_fc / num_simulations;
    else
        results.mc_collision_counts_full_covariance{trial} = [];
    end

    % Monte Carlo simulation for Sqrt QR method
    if solved_qr
        [x_hist_all_qr, u_hist_all_qr] = simulate_samples(mu_0, Sigma_0, A_sys, B_sys, G_sys, prob_qr.K, prob_qr.mu, prob_qr.v, num_nodes, num_simulations);
        collision_counts_qr = count_collision_samples(x_hist_all_qr, obstacle_centers, obstacle_radii);
        results.mc_collision_counts_qr{trial} = collision_counts_qr;
        results.mc_collision_prob_qr(trial, :, :) = collision_counts_qr / num_simulations;
    else
        results.mc_collision_counts_qr{trial} = [];
    end
end
disp("Done!")

%% 
num_samples_successfull_monte_carlo_qr = 0;
num_samples_successfull_monte_carlo_fc = 0;
for trial = 1:num_trials
    if ~isempty(results.mc_collision_counts_qr{trial}) && ...
        all(results.mc_collision_prob_qr(trial, :, :) <= state_risk, 'all')
        num_samples_successfull_monte_carlo_qr = num_samples_successfull_monte_carlo_qr + 1;
    end
    if ~isempty(results.mc_collision_counts_full_covariance{trial}) && ...
        all(results.mc_collision_prob_full_covariance(trial, :, :) <= state_risk, 'all')
        num_samples_successfull_monte_carlo_fc = num_samples_successfull_monte_carlo_fc + 1;
    end
end
fprintf('Number of successful Monte Carlo simulations for Sqrt QR: %d\n', num_samples_successfull_monte_carlo_qr);
fprintf('Number of successful Monte Carlo simulations for Full Covariance: %d\n', num_samples_successfull_monte_carlo_fc);

%% Display summary table
fprintf('\n');
fprintf('Summary table:\n');
for trial = 1:num_trials
    fprintf('  Method          | Iterations | Objective | Time (s)\n');
    fprintf('  --------------------------------------------------\n');
    fprintf('  Full Covariance | %d | %.6f | %.3f\n', results.iters_full_covariance(trial), results.objective_full_covariance(trial), results.time_full_covariance(trial));
    fprintf('  Sqrt QR         | %d | %.6f | %.3f\n', results.iters_qr(trial), results.objective_qr(trial), results.time_qr(trial));
end

%% Save all results to mat file
current_time = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
mkdir(sprintf('./data/obstacle_planning_test_%s', current_time));
save_filename = sprintf('./data/obstacle_planning_test_%s/results.mat', current_time);
fprintf('\nSaving all results to %s...\n', save_filename);
save(save_filename, 'results', 'num_trials', 'num_obstacles', 'state_risk', 'control_risk', 'num_nodes');
fprintf('Results saved successfully!\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inspect individual Monte Carlo trials
trial = 1;
prob = results.prob_full_covariance{trial};
[x_hist_all_qr, u_hist_all_qr] = simulate_samples(mu_0, Sigma_0, A_sys, B_sys, G_sys, prob.K, prob.mu, prob.v, num_nodes, num_simulations);

figure;
plot_problem(results.obstacle_centers{trial}, results.obstacle_radii{trial}, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);
plot_simulations(x_hist_all_qr);

%% Plot the results

trial = 1;
figure(Position=[0, 0, 15, 15])
sgtitle(sprintf('Trial %d', trial))

tiledlayout(2, 1);
nexttile;
title(sprintf('Full Covariance, Objective: %.3f', results.objective_full_covariance(trial)))
hold on
plot_problem(results.obstacle_centers{trial}, results.obstacle_radii{trial}, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);

if ~isempty(results.prob_full_covariance{trial}) && ~isempty(results.prob_full_covariance{trial}.mu)
	plot_solution(results.prob_full_covariance{trial}.mu, results.prob_full_covariance{trial}.P, state_risk);
end
plot(results.x_opt_deterministic{trial}(1,:), results.x_opt_deterministic{trial}(2,:), 'r.-', DisplayName='deterministic', LineWidth=1);


xlim([-1.4, 11])
ylim([-3, 2.0])

nexttile;
title(sprintf('Sqrt QR, Objective: %.3f', results.objective_qr(trial)))
hold on
plot_problem(results.obstacle_centers{trial}, results.obstacle_radii{trial}, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);

if ~isempty(results.prob_qr{trial}) && ~isempty(results.prob_qr{trial}.mu)
	plot_solution(results.prob_qr{trial}.mu, results.prob_qr{trial}.P, state_risk);
end
plot(results.x_opt_deterministic{trial}(1,:), results.x_opt_deterministic{trial}(2,:), 'r.-', DisplayName='deterministic', LineWidth=1);

% legend(legendUnq(), Location='northoutside', Orientation='horizontal', IconColumnWidth=15, FontSize=25)
xlim([-1.4, 11])
ylim([-3, 2.0])


