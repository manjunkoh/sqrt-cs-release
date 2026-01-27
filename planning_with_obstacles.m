clc; clear;
addpath(genpath('./utils'))
addpath ./src
addpath ./SCvxStar/src/
addpath ./obstacle_path_planning
figure_settings

%% Parameters
basic_parameters;

obstacle_centers = [5, 0];
obstacle_radii = 1.2;
num_obstacles = 1;

% obstacle_centers = [3, 0.5; 7.5, -1];
% obstacle_radii = [1, 1];
% num_obstacles = 2;
% 
% obstacle_centers = [2.5, -1; 4.5, 0.5; 6, -1.5; 8, -0.1];
% obstacle_radii = 0.9 * [1, 1, 1, 1];
% num_obstacles = 4;

% obstacle_centers = [7, -1.5; 1, -2; 3.5, 0.5];
% obstacle_radii = obstacle_radius * ones(1, 3);
% num_obstacles = 3;

plot_problem(obstacle_centers, obstacle_radii, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f);

% Define the chance constraints
chance_constraints_control={...
    struct('type', 'affine', 'alpha', [0; 1], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes), ...
    struct('type', 'affine', 'alpha', [0; -1], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes), ...
    struct('type', 'affine', 'alpha', [1; 0], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes), ...
    struct('type', 'affine', 'alpha', [-1; 0], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes)};

chance_constraints_state={struct('type', 'affine', 'alpha', ...
    [0; 1; 0; 0], 'beta', wall_y_pos, 'p', state_risk, 'nodes', 1:num_nodes+1)};

circular_obstacles = {};
for i = 1:num_obstacles
    circular_obstacles = [circular_obstacles, struct('center', obstacle_centers(i,:)', 'radius', obstacle_radii(i), 'p', state_risk)];
end

%% Solve the deterministic problem
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

    constraints = [constraints, x(:,1) == mu_0, x(:,num_nodes+1) == mu_f];

    constraints = [constraints, x(2,:) <= wall_y_pos];

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

%% Solve the constrained stochastic problem via iterative approach
max_iters = 100;
penalty_scalar_x = 1000;
penalty_scalar_u = 10;
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

Y_init_value = 1e-4; % This is really important!!!
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
        constraints = [constraints, Y(:,:,k) >= 0];
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
        time_full_covariance = toc;
        fprintf('Converged in %d iterations in %.3f seconds\n', iter, time_full_covariance);
        prob_fc = struct();
        prob_fc.mu = mu_opt;
        prob_fc.v = v_opt;
        prob_fc.P = P_opt;
        prob_fc.Y = Y_opt;
        prob_fc.lambda = lambda_opt;
        prob_fc.lambda_x = lambda_x_opt;

        prob.K = zeros(nu, nx, num_nodes);
        for k = 1:num_nodes
            prob_fc.K(:,:,k) = value(U(:,:,k)) / value(P(:,:,k));
        end
            
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
scp_params = SCPParams();
scp_params.k_max = 300;
scp_params.tol_opt = 1E-2;
scp_params.tol_feas = 1E-4;
scp_params.linearization = 'inexact';
% scp_params.w_inexact = 1E4;

relax_obstacle_constraints = false;

init_guess = struct();
init_guess.S = interpolate_lower_triangular(chol(Sigma_0, 'lower'), chol(Sigma_f, 'lower'), num_nodes+1, 'log-cholesky');
% init_guess.L = zeros(nu, nx, num_nodes);
init_guess.L = 1E-2 * ones(nu, nx, num_nodes);
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
    relax_obstacle_constraints=relax_obstacle_constraints);

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

fprintf('Full Covariance: Objective: %.3f, Time: %.3f seconds, Iterations: %d\n', objective_full_covariance, time_full_covariance, iters_full_covariance);
fprintf('Sqrt QR: Objective: %.3f, Time: %.3f seconds, Iterations: %d\n', objective_qr, time_qr, iters_qr);

%% Plot the covariance trajectory
% QR method
figure;
hold on
plot_problem(obstacle_centers, obstacle_radii, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);
plot_solution(prob_qr.mu, prob_qr.P, state_risk)
plot(x_opt(1,:), x_opt(2,:), 'r.-')

axis equal
xlim([-1.4, 11])
ylim([-3, 2.0])

%% Full covariance
figure;
hold on
plot_problem(obstacle_centers, obstacle_radii, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);
plot_solution(prob_fc.mu, prob_fc.P, state_risk)
plot(x_opt(1,:), x_opt(2,:), 'r.-')

axis equal
xlim([-1.4, 11])
ylim([-3, 2.0])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Evaluate the linearized constraint satisfaction for each method
% Since the linearization is a conservative approximation, violation of
% these constraints doesn't necessarily imply that the chance-constrained
% obstacle avoidance is violated. However, it should raise an alarm.

% QR method
disp("Constraint satisfaction for the QR method:")
[wall_safe_flag, obstacle_safe_flag, obstacle_safe_flags, constraint_violations_qr] ...
    = check_probabilistic_collision(prob_qr.mu, prob_qr.P, obstacle_centers, obstacle_radii, wall_y_pos, state_risk);
fprintf('Wall safe flag: %d\n', wall_safe_flag);
fprintf('Obstacle safe flag: %d\n', obstacle_safe_flag);

%% Full covariance
disp("Constraint satisfaction for the Full covariance method:")
[wall_safe_flag, obstacle_safe_flag, obstacle_safe_flags, constraint_violations_fc] ...
    = check_probabilistic_collision(prob_qr.mu, prob_qr.P, obstacle_centers, obstacle_radii, wall_y_pos, state_risk);
fprintf('Wall safe flag: %d\n', wall_safe_flag);
fprintf('Obstacle safe flag: %d\n', obstacle_safe_flag);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Monte Carlo simulation
num_simulations = 10000;

% QR method
mu = prob_qr.mu;
K = prob_qr.K;
v = prob_qr.v;

[x_hist_all, u_hist_all] = simulate_samples(mu_0, Sigma_0, A_sys, B_sys, G_sys, K, mu, v, num_nodes, num_simulations);

collision_flags_all = count_collision_samples(x_hist_all, obstacle_centers, obstacle_radii);
if all(collision_flags_all <= state_risk * num_simulations)
    fprintf('Collision probability is less than confidence level. Success.\n');
else
    fprintf('Collision probability is greater than confidence level. Failure.\n');
end

figure
plot_problem(obstacle_centers, obstacle_radii, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);
plot_simulations(x_hist_all)

%% Full covariance
mu = prob_fc.mu;
K = prob_fc.K;
v = prob_fc.v;

[x_hist_all, u_hist_all] = simulate_samples(mu_0, Sigma_0, A_sys, B_sys, G_sys, K, mu, v, num_nodes, num_simulations);

collision_flags_all = count_collision_samples(x_hist_all, obstacle_centers, obstacle_radii);
if all(collision_flags_all <= state_risk * num_simulations)
    fprintf('Collision probability is less than confidence level. Success.\n');
else
    fprintf('Collision probability is greater than confidence level. Failure.\n');
end

function control_flags_all = count_control_violations(u_hist_all, u_max)
    control_flags_all = false(size(u_hist_all));
    for i = 1:size(u_hist_all, 3)
        u_hist = u_hist_all(:,:,i);
        control_flags = u_hist > u_max | u_hist < -u_max;
        control_flags_all(:,:,i) = control_flags;
    end
    control_flags_all = sum(control_flags_all, 3);
end

control_flags_all = count_control_violations(u_hist_all, u_max);
if all(control_flags_all <= control_risk * num_simulations)
    fprintf('Control violation probability is less than confidence level. Success.\n');
else
    fprintf('Control violation probability is greater than confidence level. Failure.\n');
end

% figure
% plot_problem(obstacle_centers, obstacle_radii, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);
% plot_simulations(x_hist_all)