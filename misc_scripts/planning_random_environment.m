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

num_simulations = 10000;
%% Load generated obstacle environments
load_filename = './data/obstacle_cases_2norm.mat';
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
        
    % Create linearized hyperplane constraints
    chance_constraints_state = {};
    for k = 1:num_nodes
        for i = 1:num_obstacles
            [a, b] = hyperplane_from_circular_obstacle(obstacle_centers(i,:)', obstacle_radii(i), x_opt(pos_idx,k));
            a = [a; 0; 0];
            b = -b;
    
            chance_constraints_state = [chance_constraints_state, struct('type', 'affine', 'alpha', a, 'beta', b, 'p', state_risk, 'nodes', k)];
	    end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Solve the constrained stochastic problem via iterative approach
    max_iters = 30;
    penalty_scalar_x = 10;
    penalty_scalar_u = 100;
    penalty_increase_ratio = 2;
    max_penalty = 1e8;
    feasibility_tolerance = 1E-4;
    convergence_tolerance = 1E-3;
    
    yalmip('clear')
    mu = sdpvar(nx, num_nodes+1, 'full');
    v = sdpvar(nu, num_nodes, 'full');
    P = sdpvar(nx, nx, num_nodes+1);
    U = sdpvar(nu, nx, num_nodes, 'full');
    Y = sdpvar(nu, nu, num_nodes);
    slack_u = sdpvar(4, num_nodes, 'full'); % "virtual control"
    slack_x = sdpvar(num_obstacles, num_nodes+1, 'full');
    slack_J = sdpvar(1, num_nodes);
    
    mu_ref = x_opt;
    v_ref = u_opt;
    slack_J_ref = (u_max*0.1)^2 * ones(1, num_nodes);
    
    % Initialize tracking variables
    objective_full_covariance = NaN;
    time_full_covariance = NaN;
    
    penalty_scalar_x_prev = NaN;
    penalty_scalar_u_prev = NaN;
    objective_augmented_prev = Inf;
    
    tic;
    for iter = 1:max_iters
        constraints = [];
        objective = 0;
    
        for k = 1:num_nodes
            constraints = [constraints, mu(:,k+1) == A * mu(:,k) + B * v(:,k)];
            constraints = [constraints, P(:,:,k+1) == A * P(:,:,k) * A' + B * Y(:,:,k) * B' + A * U(:,:,k)' * B' + B * U(:,:,k) * A' + G * G'];
            constraints = [constraints, [P(:,:,k) , U(:,:,k)';
                U(:,:,k) , Y(:,:,k)] >= 0];
            constraints = [constraints, Y(:,:,k) >= 0];
        end
    
        constraints = [constraints, mu(:,1) == mu_0, mu(:,num_nodes+1) == mu_f, P(:,:,1) == Sigma_0, P(:,:,num_nodes+1) <= Sigma_f];
    
        % Control chance constraints
        for k = 1:num_nodes
            % Y_ref_k = Y_ref(:,:,k);
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
            for k = 1:num_nodes
                [a, b] = hyperplane_from_circular_obstacle(obstacle_centers(i,:)', obstacle_radii(i), x_opt(pos_idx,k));
                a = [a; 0; 0];
                b = -b;
    
                constraints = [constraints
                    z^2 * (a' * P(:,:,k) * a) <= (b - a'*mu_ref(:,k))^2 - 2 * (b - a'* mu_ref(:,k)) * a' * (mu(:,k) - mu_ref(:,k)) + slack_x(i,k)
                    b - a' * mu(:,k) >= 0
                    ];
            end
        end
    
        for k = 1:num_nodes
            constraints = [constraints
                lambda_max(Y(:,:,k)) <= slack_J_ref(k)^2 + 2 * slack_J_ref(k) * (slack_J(k) - slack_J_ref(k))
                ];
        end
    
        constraints = [constraints, slack_u(:) >= 0, slack_x(:) >= 0, slack_J >= 0];
    
        for k = 1:num_nodes
            % objective = objective + v(:,k)' * R * v(:,k) + trace(R * Y(:,:,k));
            objective = objective + norm(v(:,k)) + sqrt(chi2inv(0.99, nu)) * slack_J(k);
        end
    
        objective_augmented = objective + penalty_scalar_u * sum(slack_u, 'all') + penalty_scalar_x * sum(slack_x, 'all');
    
        fprintf("Iteration %d    ", iter)
    
        sol = optimize(constraints, objective_augmented, sdpsettings('verbose', 0));
    
    
        if sol.problem
            time_full_covariance = toc;
            iters_full_covariance = iter;
            fprintf('Infeasible at iteration %d\n', iter);
            break;
        end
    
        fprintf("Objective: %f ",  value(objective))
    
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
    
        control_constraint_values = check_control_constraint_satisfaction(v_opt, Y_opt, chance_constraints_control, num_nodes);
    
        state_constraint_values = check_state_constraint_satisfaction(mu_opt, P_opt, obstacle_centers, obstacle_radii, x_opt, pos_idx, num_nodes, num_obstacles, state_risk);
    
        fprintf("Largest cntrl violation: %f   ", max(control_constraint_values(:)));
        fprintf("Largest state violation: %f\n", max(state_constraint_values(:)));
    
        if iter > 1 ...
                && (all(slack_x_opt(:) <= feasibility_tolerance) ...
                && all(slack_u_opt(:) <= feasibility_tolerance)) ...
                && all(control_constraint_values(:) <= feasibility_tolerance) ...
                && all(state_constraint_values(:) <= feasibility_tolerance)...
                && objective_augmented_prev - objective_with_previous_penalty < convergence_tolerance
    
            time_full_covariance = toc;
            fprintf('Converged in %d iterations in %.3f seconds\n', iter, time_full_covariance);
            prob_fc = struct();
            prob_fc.mu = mu_opt;
            prob_fc.v = v_opt;
            prob_fc.P = P_opt;
            prob_fc.Y = Y_opt;
            prob_fc.P_u = Y_opt;
            prob_fc.lambda_u = slack_u_opt;
            prob_fc.lambda_x = slack_x_opt;
    
            prob_fc.dv99 = 0;
            for k = 1:num_nodes
                prob_fc.dv99 = prob_fc.dv99 + norm(v_opt) + sqrt(chi2inv(0.99, nu)) * sqrt(lambda_max(Y_opt(:,:,k)));
            end
    
            prob_fc.K = zeros(nu, nx, num_nodes);
            for k = 1:num_nodes
                prob_fc.K(:,:,k) = value(U(:,:,k)) / value(P(:,:,k));
            end
    
            break;
        end
    
        v_ref = v_opt;
        mu_ref = mu_opt;
    
        penalty_scalar_x_prev = penalty_scalar_x;
        penalty_scalar_u_prev = penalty_scalar_u;
        penalty_scalar_x = min(max_penalty, penalty_scalar_x * penalty_increase_ratio);
        penalty_scalar_u = min(max_penalty, penalty_scalar_u * penalty_increase_ratio);
        objective_augmented_prev = objective_augmented;
    
    
        if iter == max_iters
            time_full_covariance = toc;
            fprintf('Reached maximum iterations\n')
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
        objective_type='DV99', ...
        mu_0=mu_0, mu_f=mu_f, ...
        P_0=Sigma_0, P_f=Sigma_f, Q=Q, R=R, ...
        chance_constraints_state=chance_constraints_state, ...
        chance_constraints_control=chance_constraints_control);

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
trial = 3;
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


