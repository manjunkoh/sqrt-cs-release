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
    struct('type', 'affine', 'alpha', [1; 0], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes), ...
    struct('type', 'affine', 'alpha', [-1; 0], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes), ...
    struct('type', 'affine', 'alpha', [0; 1], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes), ...
    struct('type', 'affine', 'alpha', [0; -1], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes)};

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
slack_u = sdpvar(1, num_nodes, 'full');

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
                a' * x(pos_idx,k) + b <= slack_u(k)
            ];
        end
    end

    constraints = [constraints, slack_u >= 0];

    for k = 1:num_nodes
        objective = objective + norm(u(:,k)) + penalty_scalar_obstacle * slack_u(k);
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

%%
figure;
hold on
plot_problem(obstacle_centers, obstacle_radii, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);
plot(x_opt(1,:), x_opt(2,:), 'r.-')

%% Generate hyperplane constraints for each node&obstacle based on the solution from the deterministic problem
% chance_constraints_state={struct('type', 'affine', 'alpha', ...
    % [0; 1; 0; 0], 'beta', wall_y_pos, 'p', state_risk, 'nodes', 1:num_nodes+1)};
chance_constraints_state = {};

for k = 1:num_nodes
    for i = 1:num_obstacles
        [a, b] = hyperplane_from_circular_obstacle(obstacle_centers(i,:)', obstacle_radii(i), x_opt(pos_idx,k));
        a = [a; 0; 0];
        b = -b;

        chance_constraints_state = [chance_constraints_state, struct('type', 'affine', 'alpha', a, 'beta', b, 'p', state_risk, 'nodes', k)];
	end
end


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
init_guess.L = zeros(nu, nx, num_nodes);
init_guess.mu = x_opt;
init_guess.v = u_opt;

prob_qr = SqrtQRCovarianceSteering(init_guess, ...
	N=num_nodes, ...
	A_sys=A_sys, B_sys=B_sys, G_sys=G_sys, ...
    objective_type = 'DV99', ...
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
else
	time_qr = NaN;
	iters_qr = NaN;
	objective_qr = NaN;
	fprintf('SqrtQRCovarianceSteering did not solve successfully\n');
end

%%
figure;
hold on
plot_problem(obstacle_centers, obstacle_radii, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);
plot_solution(prob_qr.mu, prob_qr.P, state_risk)
plot(x_opt(1,:), x_opt(2,:), 'r.-')
% plot(x_opt(1,1:15), x_opt(2,1:15), 'b.-')
% plot(prob_qr.mu(1,1:15), prob_qr.mu(2,1:15), 'g.-')
axis equal
% xlim([-1.4, 11])
% ylim([-3, 2.0])

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
		&& all(state_constraint_values(:) <= feasibility_tolerance)
        % && objective_augmented_prev - objective_with_previous_penalty < convergence_tolerance
        
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

        prob.K = zeros(nu, nx, num_nodes);
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

%% Full covariance
figure;
hold on
plot_problem(obstacle_centers, obstacle_radii, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);
plot_solution(prob_fc.mu, prob_fc.P, state_risk)
plot(x_opt(1,:), x_opt(2,:), 'r.-')

axis equal
xlim([-1.4, 11])
ylim([-3, 2.0])

%% Monte Carlo simulation
rng(1)
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

control_flags_all = count_control_violations(u_hist_all, u_max);
if all(control_flags_all <= control_risk * num_simulations)
    fprintf('Control violation probability is less than confidence level. Success.\n');
else
    fprintf('Control violation probability is greater than confidence level. Failure.\n');
end

% figure
% plot_problem(obstacle_centers, obstacle_radii, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);
% plot_simulations(x_hist_all)

%%
figure;
tiledlayout(nu, 1)
for i = 1:nu
    nexttile
	hold on
	for sample_idx = 1:num_simulations
		plot(u_hist_all(i,:,sample_idx));
	end
end

%% Full covariance
num_simulations = 10000;
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

control_flags_all = count_control_violations(u_hist_all, u_max);
if all(control_flags_all <= control_risk * num_simulations)
    fprintf('Control violation probability is less than confidence level. Success.\n');
else
    fprintf('Control violation probability is greater than confidence level. Failure.\n');
end


function control_flags_all = count_control_violations(u_hist_all, u_max)
    control_flags_all = false(4, size(u_hist_all, 2), size(u_hist_all, 3));
    for constraint_idx = 1:4
		switch constraint_idx
			case 1
				control_flags = u_hist_all(1,:,:) > u_max;
			case 2
				control_flags = u_hist_all(1,:,:) < -u_max;
			case 3
				control_flags = u_hist_all(2,:,:) > u_max;
			case 4
				control_flags = u_hist_all(2,:,:) < -u_max;
		end
        control_flags_all(constraint_idx,:,:) = control_flags;
    end
    control_flags_all = sum(control_flags_all, 3);
end

%%
figure;
tiledlayout(nu, 1)
for i = 1:nu
    nexttile
	hold on
	for sample_idx = 1:num_simulations
		plot(u_hist_all(i,:,sample_idx));
	end
end

%% Check constraint satisfaction

control_constraint_values = check_control_constraint_satisfaction(v_opt, Y_opt, chance_constraints_control, num_nodes)


function constraint_values = check_control_constraint_satisfaction(v, Y, chance_constraints_control, num_nodes)
    constraint_values = NaN(4, num_nodes);
    for k = 1:num_nodes
	    for i = 1:length(chance_constraints_control)
		    a = chance_constraints_control{i}.alpha;
		    b = chance_constraints_control{i}.beta;
		    p = chance_constraints_control{i}.p;
		    z = norminv(1 - p);
		    constraint_values(i,k) = [
			    z^2 * (a' * Y(:,:,k) * a) - (b - a'*v(:,k))^2
			    ];
	    end
    end
	constraint_values = max(constraint_values, 0);
end


function constraint_values = check_state_constraint_satisfaction(mu, P, obstacle_centers, obstacle_radii, x_opt, pos_idx, num_nodes, num_obstacles, state_risk, mu_ref)
    constraint_values = NaN(num_obstacles, num_nodes);
	constraint_values_linearized = NaN(num_obstacles, num_nodes);
	z = norminv(1 - state_risk);
    for k = 1:num_nodes+1
        for i = 1:num_obstacles
            [a, b] = hyperplane_from_circular_obstacle(obstacle_centers(i,:)', obstacle_radii(i), x_opt(pos_idx,k));
            a = [a; 0; 0];
            b = -b;
            constraint_values(i,k) = [
                z^2 * (a' * P(:,:,k) * a) - (b - a'*mu(:,k))^2
            ];
			% constraint_values_linearized(i,k) = [
				% z^2 * (a' * P(:,:,k) * a) - (b - a'*mu_ref(pos_idx,k))^2 - 
        end
    end
	constraint_values = max(constraint_values, 0);
end

%% Plot the state components' sigmas
formulations = {prob_fc, prob_qr};
formulation_names = {'FullCov', 'SqrtQR'};
figure;
tiledlayout(nx, 1);
for i = 1:nx
	nexttile;
    hold on
	for formulation_idx = 1:length(formulations)
		formulation = formulations{formulation_idx};
		one_sigma_plus = formulation.mu(i,:) + sqrt(squeeze(formulation.P(i,i,:)))';
		% one_sigma_minus = formulation.mu(i,:) - sqrt(squeeze(formulation.P(i,i,:))');
		plot(0:num_nodes, one_sigma_plus, DisplayName=sprintf('%s', formulation_names{formulation_idx}));
		% plot(0:num_nodes, one_sigma_minus, 'b-', DisplayName=sprintf('%s', formulation_names{formulation_idx}));
	end
	legend(legendUnq(), Location='northoutside', Orientation='horizontal', IconColumnWidth=15, FontSize=25)
end

%% Plot the control components' sigmas
figure;
tiledlayout(nu, 1);
for i = 1:nu
	nexttile;
	hold on
	for formulation_idx = 1:length(formulations)
		formulation = formulations{formulation_idx};
		one_sigma_plus = formulation.v(i,:) + sqrt(squeeze(formulation.P_u(i,i,:)))';
		% one_sigma_minus = formulation.v(i,:) - sqrt(squeeze(formulation.P_u(i,i,:))');
		stairs(0:num_nodes-1, one_sigma_plus, DisplayName=sprintf('%s', formulation_names{formulation_idx}));
		% plot(0:num_nodes, one_sigma_minus, 'b-', DisplayName=sprintf('%s', formulation_names{formulation_idx}));
	end
	legend(legendUnq(), Location='northoutside', Orientation='horizontal', IconColumnWidth=15, FontSize=25)
end

