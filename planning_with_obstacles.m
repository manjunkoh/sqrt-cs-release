clc; clear;
addpath(genpath('./utils'))
addpath ./src
addpath ./SCvxStar/src/
figure_settings

function plot_obstacles(centers, radii)
    hold on
    for i = 1:size(centers, 1)
        filled_circle(centers(i, :), radii(i));
        text(centers(i, 1), centers(i, 2), sprintf('%d', i));
    end
end

function filled_circle(center, radius)
    rectangle('Position', [center(1)-radius, center(2)-radius, 2*radius, 2*radius], 'Curvature', [1, 1], 'FaceColor', 'k', 'FaceAlpha', 0.5);
end

function plot_wall(wall_y_pos)
    yline(wall_y_pos, 'k');
end

function plot_problem(obstacle_center, obstacle_radii, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, options)
    arguments
        obstacle_center
        obstacle_radii
        wall_y_pos
        mu_0
        Sigma_0
        mu_f
        Sigma_f
        options.fig = []
    end

    if isempty(options.fig)
        figure(Position=[0, 0, 20, 10]);
    else
        figure(options.fig);
    end
    axis equal
    plot_obstacles(obstacle_center, obstacle_radii);
    plot_wall(wall_y_pos);
    plot3sigmaEllipse(mu_0, Sigma_0, Color='#D55E00', DisplayName='Start')
    plot3sigmaEllipse(mu_f, Sigma_f, Color='#D55E00', LineStyle=":", DisplayName='Goal')
    xlabel('$x$')
    ylabel('$y$')
end

function flags = is_collision_free_per_node_per_obstacle(x_hist, obstacle_centers, obstacle_radii)
    flags = true(size(x_hist, 2), size(obstacle_centers, 1));
    for k = 1:size(x_hist, 2)
        for i = 1:size(obstacle_centers, 1)
            if norm(x_hist(1:2,k) - obstacle_centers(i,:)') < obstacle_radii(i)
                flags(k,i) = false;
            end
        end
    end
end

function yn = is_collision_free(x_opt, obstacle_centers, obstacle_radii)
    yn = false;
    for k = 1:size(x_opt, 2)
        for i = 1:size(obstacle_centers, 1)
            if norm(x_opt(1:2,k) - obstacle_centers(i,:)') < obstacle_radii(i)
                return
            end
        end
    end
    yn = true;
end

function is_valid = is_valid_environment(obstacle_centers, obstacle_radii, mu_0, mu_f, Sigma_0, Sigma_f)
    % if any of the obstacles are too close to the start or goal, return false
    is_valid = false;
    for i = 1:size(obstacle_centers, 1)
        if norm(mu_0(1:2) - obstacle_centers(i,:)') < obstacle_radii(i) + 3 * sqrt(Sigma_0(1,1))
            return
        end
        if norm(mu_f(1:2) - obstacle_centers(i,:)') < obstacle_radii(i) + 3 * sqrt(Sigma_f(1,1))
            return
        end
    end
    is_valid = true;
end

function [obstacle_centers, obstacle_radii] = generate_random_environment(num_obstacles, obstacle_radius, wall_y_pos, mu_0, mu_f, Sigma_0, Sigma_f)
    is_valid = false;
    while ~is_valid
        obstacle_centers = [rand(num_obstacles, 1) * 10, rand(num_obstacles, 1) * (wall_y_pos - (-3)) + (-3)];
        obstacle_radii = obstacle_radius * ones(1, num_obstacles);
        is_valid = is_valid_environment(obstacle_centers, obstacle_radii, mu_0, mu_f, Sigma_0, Sigma_f);
    end
end

%% Parameters
wall_y_pos = 1.5;
% obstacle_centers = [3, 0.5; 7.5, -1];
% obstacle_radii = [1, 1];

% obstacle_centers = [2.5, -1; 4.5, 0.5; 6, -1.5; 8, -0.1];
% obstacle_radius = 1;

% num_obstacles = size(obstacle_centers, 1);
pos_idx = 1:2;

num_nodes = 50;
t_f = 30;
dt = t_f / num_nodes;
nx = 4;
nu = 2;

mu_0 = [0; 0; 0; 0];
mu_f = [10; 0; 0; 0];

Sigma_0 = diag([0.1, 0.1, 0.01, 0.01]);
% Sigma_f = Sigma_0;
Sigma_f = diag([0.04, 0.04, 0.01, 0.01]);

Q = 0.001 * eye(nx);
R = 0.01 * eye(nu);

A = [eye(2), dt*eye(2); zeros(2, 2), eye(2)];
B = [0.5*dt^2*eye(2); dt*eye(2)];
% G = 0.01 * [zeros(2, 2); eye(2)];
% G = [0.01 * eye(2); 0.01 * eye(2)];
q = 0.005;
G = sqrt(q * dt) * eye(4);

A_sys = repmat(A, [1,1,num_nodes]);
B_sys = repmat(B, [1,1,num_nodes]);
G_sys = repmat(G, [1,1,num_nodes]);

control_risk = 0.005;
u_max = 0.15;
state_risk = 0.005;

% Define the chance constraints
chance_constraints_control={...
    struct('type', 'affine', 'alpha', [0; 1], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes), ...
    struct('type', 'affine', 'alpha', [0; -1], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes), ...
    struct('type', 'affine', 'alpha', [1; 0], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes), ...
    struct('type', 'affine', 'alpha', [-1; 0], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes)};

chance_constraints_state={struct('type', 'affine', 'alpha', ...
    [0; 1; 0; 0], 'beta', wall_y_pos, 'p', state_risk, 'nodes', 1:num_nodes+1)};


num_trials = 10;
results = struct();
results.iters_full_covariance = zeros(num_trials, 1);
results.objective_full_covariance = zeros(num_trials, 1);
results.time_full_covariance = zeros(num_trials, 1);
results.prob_full_covariance = cell(num_trials, 1);
results.iters_qr = zeros(num_trials, 1);
results.objective_qr = zeros(num_trials, 1);
results.time_qr = zeros(num_trials, 1);
results.prob_qr = cell(num_trials, 1);

num_obstacles = 3;
obstacle_radius = 1;

rng(0)

for trial = 1:num_trials
    fprintf('Trial %d/%d ... ', trial, num_trials);
    [obstacle_centers, obstacle_radii] = generate_random_environment(num_obstacles, obstacle_radius, wall_y_pos, mu_0, mu_f, Sigma_0, Sigma_f);

    circular_obstacles = {};
    for i = 1:num_obstacles
        circular_obstacles = [circular_obstacles, struct('center', obstacle_centers(i,:)', 'radius', obstacle_radii(i), 'p', state_risk)];
    end

    results.obstacle_centers{trial} = obstacle_centers;
    results.obstacle_radii{trial} = obstacle_radii;
    % plot_scenario = @()plot_problem(obstacle_centers, obstacle_radii, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f);
    % plot_scenario()

    %% Solve the deterministic problem
    flag_deterministic = false;
    penalty_scalar_obstacle = 100;
    max_iters = 50;
    
    yalmip('clear')
    x = sdpvar(nx, num_nodes+1, 'full');
    u = sdpvar(nu, num_nodes, 'full');
    lambda = sdpvar(1, num_nodes, 'full');
    
    x_ref = linspace_vec(mu_0, mu_f, num_nodes+1);
    
    for iter = 1:max_iters
    
        fprintf("Iter %d\n", iter)
    
        constraints = [];
        objective = 0;
    
        for k = 1:num_nodes
            constraints = [constraints, x(:,k+1) == A * x(:,k) + B * u(:,k)];
            constraints = [constraints, u(:,k) <= u_max];
            constraints = [constraints, u(:,k) >= -u_max];
        end
        % if iter == 1
            % constraints = [constraints, x(2, ceil(num_nodes* 0.5)) >= 0.1]; % symmetry-breaking constraint
        % end
        constraints = [constraints, x(:,1) == mu_0, x(:,num_nodes+1) == mu_f];
    
        for k = 1:num_nodes
            for i = 1:size(obstacle_centers, 1)
                a = - (x_ref(pos_idx,k) - obstacle_centers(i,:)');
                b = - 0.5 * norm(x_ref(pos_idx,k) - obstacle_centers(i,:)')^2  + 0.5 * obstacle_radii(i)^2 - a' * x_ref(pos_idx,k);
                constraints = [constraints
                    a' * x(pos_idx,k) + b <= lambda(k)
                ];
            end
        end
    
        constraints = [constraints, lambda >= 0];
    
        for k = 1:num_nodes
            objective = objective + x(:,k)' * Q * x(:,k) + u(:,k)' * R * u(:,k) + penalty_scalar_obstacle * lambda(k);
        end
    
        sol = optimize(constraints, objective, sdpsettings('verbose', 0));
    
        if sol.problem
            disp('Constrained problem infeasible');
            break;
        end
    
        x_opt = value(x);
        u_opt = value(u);
    
        if norm(x_opt - x_ref) < 1e-3 && is_collision_free(x_opt, obstacle_centers, obstacle_radii)
            flag_deterministic = true;
            fprintf('Converged in %d iterations\n', iter);
            break;
        end
    
        x_ref = x_opt;
    
        if iter == max_iters
            fprintf('Reached max iters.\n')
        end
    
    end

    if ~flag_deterministic
        fprintf('Deterministic problem unsolved. Skipping stochastic problem.\n');
        results.iters_full_covariance(trial) = NaN;
        results.objective_full_covariance(trial) = NaN;
        results.time_full_covariance(trial) = NaN;
        results.iters_qr(trial) = NaN;
        results.objective_qr(trial) = NaN;
        results.time_qr(trial) = NaN;
        results.prob_full_covariance{trial} = [];
        results.prob_qr{trial} = [];
        continue;
    end

    %% Solve the constrained stochastic problem via iterative approach
    max_iters = 100;
    penalty_scalar_x = 10000;
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

    Y_init_value = 1e-3;
    mu_ref = x_opt;
    v_ref = u_opt;
    % mu_ref = linspace_vec(mu_0, mu_f, num_nodes+1);
    % v_ref = zeros(nu, num_nodes);

    P_ref = repmat(Sigma_0, [1,1,num_nodes+1]);
    Y_ref = repmat(Y_init_value * eye(nu), [1,1,num_nodes]); % This is really important!!!

    % Initialize tracking variables
    iters_full_covariance = NaN;
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
            iters_full_covariance = iter;
            fprintf('Infeasible at iteration %d\n', iter);
            break;
        end

        fprintf("Objective: %f\n",  value(objective))

        mu_opt = value(mu);
        v_opt = value(v);
        P_opt = value(P);
        Y_opt = value(Y);
        lambda_opt = value(lambda);
        lambda_x_opt = value(lambda_x);
        
        objective_full_covariance = value(objective);

        if all(vecnorm(mu_opt - mu_ref, Inf) < convergence_tolerance) 
            % && is_collision_free(mu_opt, obstacle_centers, obstacle_radius)
            time_full_covariance = toc;
            fprintf('Converged in %d iterations in %.3f seconds\n', iter, time_full_covariance);
            prob_fc = struct();
            prob_fc.mu = mu_opt;
            prob_fc.v = v_opt;
            prob_fc.P = P_opt;
            prob_fc.Y = Y_opt;
            prob_fc.lambda = lambda_opt;
            prob_fc.lambda_x = lambda_x_opt;
            break;
        end

        mu_ref = mu_opt;
        P_ref = P_opt;

        if iter == max_iters
            time_full_covariance = toc;
            fprintf('Reached maximum iterations\n')
        end
    end

    iters_full_covariance = iter;

    %% Solve with SQRT QR method

    init_guess = struct();
    init_guess.S = interpolate_lower_triangular(chol(Sigma_0, 'lower'), chol(Sigma_f, 'lower'), num_nodes+1, 'log-cholesky');
    init_guess.L = zeros(nu, nx, num_nodes);
    init_guess.mu = x_opt;
    init_guess.v = u_opt;
    % init_guess.mu = linspace_vec(mu_0, mu_f, num_nodes+1);
    % init_guess.v = zeros(nu, num_nodes);


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

    scp_params = SCPParams();
    scp_params.k_max = 100;
    scp_params.tol_opt = 1E-2;
    scp_params.tol_feas = 1E-4;
    % scp_params.w_init = 1000;
    % scp_params.w_inexact = 10;
    % scp_params.r_init = 0.5;
    scp_params.linearization = 'inexact';

    flag_solved_qr = prob_qr.solve(scp_params=scp_params);

    if flag_solved_qr
        prob_qr.postprocess();
        time_qr = seconds(prob_qr.scp.report.time);
        iters_qr = prob_qr.scp.report.iters;
        objective_qr = prob_qr.objective(prob_qr.sol);
        fprintf('SqrtQRCovarianceSteering solved successfully in %.3f seconds\n', time_qr);
        fprintf('  Number of iterations: %d\n', iters_qr);
        fprintf('  Objective value: %.6f\n', objective_qr);
    else
        time_qr = NaN;
        iters_qr = NaN;
        objective_qr = NaN;
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

%% Display summary table
fprintf('\n');
fprintf('Summary table:\n');
for trial = 1:num_trials
    fprintf('  Method          | Iterations | Objective | Time (s)\n');
    fprintf('  --------------------------------------------------\n');
    fprintf('  Full Covariance | %d | %.6f | %.3f\n', results.iters_full_covariance(trial), results.objective_full_covariance(trial), results.time_full_covariance(trial));
    fprintf('  Sqrt QR         | %d | %.6f | %.3f\n', results.iters_qr(trial), results.objective_qr(trial), results.time_qr(trial));
end

%% Plot the results
current_time = datetime('now', 'Format', 'yyyy-MM-dd_HH-mm-ss');
mkdir(sprintf('./data/obstacle_planning_test_%s', current_time));
save_filename = sprintf('./data/obstacle_planning_test_%s/results.mat', current_time);
% save(save_filename, 'results');

for trial = 1:num_trials
    figure(Position=[0, 0, 15, 15])
    sgtitle(sprintf('Trial %d', trial))

    tiledlayout(2, 1);
    nexttile;
    title(sprintf('Full Covariance, Objective: %.3f', results.objective_full_covariance(trial)))
    hold on
    plot_problem(results.obstacle_centers{trial}, results.obstacle_radii{trial}, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);

    if ~isempty(results.prob_full_covariance{trial}) && ~isempty(results.prob_full_covariance{trial}.mu)
        plot_solution(results.prob_full_covariance{trial}.mu, results.prob_full_covariance{trial}.P, state_risk);
        plot(x_opt(1,:), x_opt(2,:), 'r.-', DisplayName='deterministic', LineWidth=1);
    end

    xlim([-1.4, 11])
    ylim([-3, 2.0])

    nexttile;
    title(sprintf('Sqrt QR, Objective: %.3f', results.objective_qr(trial)))
    hold on
    plot_problem(results.obstacle_centers{trial}, results.obstacle_radii{trial}, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);

    if ~isempty(results.prob_qr{trial}) && ~isempty(results.prob_qr{trial}.mu)
        plot_solution(results.prob_qr{trial}.mu, results.prob_qr{trial}.P, state_risk);
        plot(x_opt(1,:), x_opt(2,:), 'r.-', DisplayName='deterministic', LineWidth=1);
    end
    % legend(legendUnq(), Location='northoutside', Orientation='horizontal', IconColumnWidth=15, FontSize=25)
    xlim([-1.4, 11])
    ylim([-3, 2.0])

    save_filename = sprintf('./data/obstacle_planning_test_%s/trial_%d.pdf', current_time, trial);
    exportgraphics(gcf, save_filename, ContentType='vector');

end

function plot_solution(mu, P, state_risk)
    for k = 1:size(mu, 2)
        fill3sigmaEllipse(mu(:,k), P(:,:,k), '', FaceColor='#0082B2', FaceAlpha=0.5, EdgeColor='none', DisplayName="$3 \sigma$ ellipse");
    end
    plot(mu(1,:), mu(2,:), 'k.-', DisplayName='mean', LineWidth=1);
end

%% Monte Carlo simulation
num_trial = 3;
num_simulations = 10000;

mu = results.prob_qr{num_trial}.mu;
K = results.prob_qr{num_trial}.K;
v = results.prob_qr{num_trial}.v;
obstacle_centers = results.obstacle_centers{num_trial};
obstacle_radii = results.obstacle_radii{num_trial};

[x_hist_all, u_hist_all] = simulate_samples(mu_0, Sigma_0, A_sys, B_sys, G_sys, K, mu, v, num_nodes, num_simulations);

collision_flags_all = count_collision_samples(x_hist_all, obstacle_centers, obstacle_radii);

figure
plot_problem(results.obstacle_centers{num_trial}, results.obstacle_radii{num_trial}, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);
plot_simulations(x_hist_all)

function [x_hist_all, u_hist_all] = simulate_samples(mu_0, Sigma_0, A_sys, B_sys, G_sys, K, mu, v, num_nodes, num_simulations)
    disp('Simulating samples...');
    x_hist_all = zeros(size(mu, 1), num_nodes+1, num_simulations);
    u_hist_all = zeros(size(v, 1), num_nodes, num_simulations);
    for i = 1:num_simulations
        [x_hist, u_hist] = simulate_sample(mu_0, Sigma_0, A_sys, B_sys, G_sys, K, mu, v, num_nodes);
        x_hist_all(:,:,i) = x_hist;
        u_hist_all(:,:,i) = u_hist;
    end
    disp('Done simulating samples.');
end

function [x_hist, u_hist] = simulate_sample(mu_0, Sigma_0, A_sys, B_sys, G_sys, K, mu, v, num_nodes)
    x_hist = zeros(size(mu_0, 1), num_nodes+1);
    u_hist = zeros(size(v, 1), num_nodes);
    x_k = mu_0 + chol(Sigma_0, 'lower') * randn(size(mu_0));
    x_hist(:,1) = x_k;
    for k = 1:num_nodes
        u_k = v(:,k) + K(:,:,k) * (x_k - mu(:,k));
        x_k = A_sys(:,:,k) * x_k + B_sys(:,:,k) * u_k + G_sys(:,:,k) * randn(size(G_sys, 2), 1);

        x_hist(:,k+1) = x_k;
        u_hist(:,k) = u_k;
    end
end

function plot_simulations(x_hist_all)
    for i = 1:size(x_hist_all, 3)
        x_hist = x_hist_all(:,:,i);
        plot(x_hist(1,:), x_hist(2,:), Color='g', LineWidth=0.5);
    end
end

function [collision_flags_all] = count_collision_samples(x_hist_all, obstacle_centers, obstacle_radii)
    collision_flags_all = false(size(x_hist_all, 2), size(obstacle_centers, 1), size(x_hist_all, 3));
    for i = 1:size(x_hist_all, 3)
        x_hist = x_hist_all(:,:,i);
        collision_flags = ~is_collision_free_per_node_per_obstacle(x_hist, obstacle_centers, obstacle_radii);
        collision_flags_all(:,:,i) = collision_flags;
    end

    collision_flags_all = sum(collision_flags_all, 3);
end


%% Tests

num_trial = 1;
[wall_safe_flag, obstacle_safe_flag] = check_probabilistic_collision(results.prob_full_covariance{num_trial}.mu, results.prob_full_covariance{num_trial}.P, results.obstacle_centers{num_trial}, results.obstacle_radii{num_trial}, wall_y_pos, state_risk);
fprintf('Wall safe flag: %d\n', wall_safe_flag);
fprintf('Obstacle safe flag: %d\n', obstacle_safe_flag);


%%
num_trial = 1;
[wall_safe_flag, obstacle_safe_flag, obstacle_safe_flags] = check_probabilistic_collision(results.prob_qr{num_trial}.mu, results.prob_qr{num_trial}.P, results.obstacle_centers{num_trial}, results.obstacle_radii{num_trial}, wall_y_pos, state_risk);
fprintf('Wall safe flag: %d\n', wall_safe_flag);
fprintf('Obstacle safe flag: %d\n', obstacle_safe_flag);

function [wall_safe_flag, obstacle_safe_flag, obstacle_safe_flags] = check_probabilistic_collision(mu, P, obstacle_centers, obstacle_radii, wall_y_pos, state_risk)

    if isempty(mu) || isempty(P)
        wall_safe_flag = NaN;
        obstacle_safe_flag = NaN;
        return;
    end
    wall_safe_flag = true;
    obstacle_safe_flags = true(size(obstacle_centers, 1), size(mu, 2));
    constraint_violations = zeros(size(obstacle_centers, 1), size(mu, 2));
    
    pos_idx = 1:2;
    z = norminv(1 - state_risk);
    for k = 1:size(mu, 2)
        for i = 1:size(obstacle_centers, 1)
            % consider a hyperplane that passes through the obstacle center and is perpendicular to the x-axis
            [a, b] = hyperplane_from_circular_obstacle(obstacle_centers(i,:)', obstacle_radii(i), mu(pos_idx,k));
            constraint_violations(i,k) = a' * mu(pos_idx,k) + b + z * sqrt(a' * P(pos_idx,pos_idx,k) * a);
            obstacle_safe_flags(i,k) = constraint_violations(i,k) <= 0;
        end

        wall_safe_flag = wall_safe_flag && (mu(2,k) - wall_y_pos + z * sqrt(P(2,2,k)) <= 0);
    end

    constraint_violations = max(0, constraint_violations);

    obstacle_safe_flag = all(obstacle_safe_flags, 'all');


    figure(Position=[0, 0, 15, 15]);
    hold on

    plot_obstacles(obstacle_centers, obstacle_radii);
    plot_wall(wall_y_pos);

    axis equal

    % xlim([-1.4, 11])
    % ylim([-3, 2.0])

    for k = 1:size(mu, 2)
        % for i = 1:size(obstacle_centers, 1)
        for i = 2
            % if constraint_violations(i,k) > 0
                if constraint_violations(i,k) > 0
                    FaceColor = 'r';
                else
                    FaceColor = '#0082B2';
                end

                fill_confidence_ellipse(mu(pos_idx,k), P(pos_idx,pos_idx,k), 1-state_risk, '', FaceColor=FaceColor, FaceAlpha=0.5, EdgeColor='none', DisplayName="$3 \sigma$ ellipse");

                [a, b] = hyperplane_from_circular_obstacle(obstacle_centers(i,:)', obstacle_radii(i), mu(pos_idx,k));
                constraint_violation = a' * mu(pos_idx,k) + b + z * sqrt(a' * P(pos_idx, pos_idx, k) * a);
                % plot_hyperplane(a, b, [3, 7], 'k-', LineWidth=1);
                % line([obstacle_centers(i,1), mu(1,k)], [obstacle_centers(i,2), mu(2,k)], 'Color', 'r', 'LineWidth', 1);
                drawnow
                keyboard
            % end

        end
    end

end

            

%%
num_trial = 1;
figure
hold on
plot_problem(results.obstacle_centers{num_trial}, results.obstacle_radii{num_trial}, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);
% plot_solution(results.prob_qr{num_trial}.mu, results.prob_qr{num_trial}.P, state_risk);
plot(results.prob_qr{num_trial}.mu(1,:), results.prob_qr{num_trial}.mu(2,:), 'r.-')
plot(results.prob_qr{num_trial}.scp.this_iter.ref_vars.mu(1,:), results.prob_qr{num_trial}.scp.this_iter.ref_vars.mu(2,:), 'b.-')

%%
function plot_hyperplane(a, b, xlim, varargin)
    x = linspace(xlim(1), xlim(2), 100);
    y = (-a(1) * x - b) / a(2);
    plot(x, y, varargin{:});
end

p = 0.99;
mu = [0, 1]';
a = [0, 1]';
b = -4;
P = [1, 0; 0, 0.09];

% y = mu - a / norm(a) * norminv(p) * sqrt(a' * P * a);

constraintLHS = a' * mu + b + norminv(p) * sqrt(a' * P * a) % want this to be <= 0

figure
hold on
% plot1sigmaEllipse(mu, P)
plot_confidence_ellipse(mu, P, p)
plot_hyperplane(a, b, [-2, 3])
% plot(y(1), y(2), 'g.')

axis equal