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

plot_problem(obstacle_centers, obstacle_radii, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f);

% Define the chance constraints
chance_constraints_control={...
    struct('type', 'affine', 'alpha', [1; 0], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes), ...
    struct('type', 'affine', 'alpha', [-1; 0], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes), ...
    struct('type', 'affine', 'alpha', [0; 1], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes), ...
    struct('type', 'affine', 'alpha', [0; -1], 'beta', u_max, 'p', control_risk, 'nodes', 1:num_nodes)};

%% Solve the deterministic problem
slack_penalty_obstacle = 100;
max_iters = 50;

yalmip('clear')
x = sdpvar(nx, num_nodes+1, 'full');
u = sdpvar(nu, num_nodes, 'full');
slack = sdpvar(1, num_nodes, 'full');

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
                a' * x(pos_idx,k) + b <= slack(k)
            ];
        end
    end

    constraints = [constraints, slack >= 0];

    for k = 1:num_nodes
        objective = objective + norm(u(:,k)) + slack_penalty_obstacle * slack(k);
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

%% Solve with block method
prob_block = BlockCholeskySteering(...
    A=A_sys, B=B_sys, G=G_sys, ...
    P_0=Sigma_0, P_f=Sigma_f, ...
    objective_type='DV99', ...
    mu_0=mu_0, mu_f=mu_f, ...
    chance_constraints_state=chance_constraints_state, ...
    chance_constraints_control=chance_constraints_control, ...
    N=num_nodes);

t_block = tic;
diag_block = prob_block.solve(sdpsettings('verbose', 1));

%% Solve with SQRT QR method
scp_params = SCPParams();
scp_params.k_max = 300;
scp_params.tol_opt = 1E-2;
scp_params.tol_feas = 1E-4;

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
axis equal
% xlim([-1.4, 11])
% ylim([-3, 2.0])

%% Solve the constrained stochastic problem via convex-concave procedure
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
% xlim([-1.4, 11])
% ylim([-3, 2.0])


% figure
% plot_problem(obstacle_centers, obstacle_radii, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);
% plot_simulations(x_hist_all)


% %% Plot the state components' sigmas
% formulations = {prob_fc, prob_qr};
% formulation_names = {'FullCov', 'SqrtQR'};
% figure;
% tiledlayout(nx, 1);
% for i = 1:nx
% 	nexttile;
%     hold on
% 	for formulation_idx = 1:length(formulations)
% 		formulation = formulations{formulation_idx};
% 		one_sigma_plus = formulation.mu(i,:) + sqrt(squeeze(formulation.P(i,i,:)))';
% 		% one_sigma_minus = formulation.mu(i,:) - sqrt(squeeze(formulation.P(i,i,:))');
% 		plot(0:num_nodes, one_sigma_plus, DisplayName=sprintf('%s', formulation_names{formulation_idx}));
% 		% plot(0:num_nodes, one_sigma_minus, 'b-', DisplayName=sprintf('%s', formulation_names{formulation_idx}));
% 	end
% 	legend(legendUnq(), Location='northoutside', Orientation='horizontal', IconColumnWidth=15, FontSize=25)
% end

% %% Plot the control components' sigmas
% figure;
% tiledlayout(nu, 1);
% for i = 1:nu
% 	nexttile;
% 	hold on
% 	for formulation_idx = 1:length(formulations)
% 		formulation = formulations{formulation_idx};
% 		one_sigma_plus = formulation.v(i,:) + sqrt(squeeze(formulation.P_u(i,i,:)))';
% 		% one_sigma_minus = formulation.v(i,:) - sqrt(squeeze(formulation.P_u(i,i,:))');
% 		stairs(0:num_nodes-1, one_sigma_plus, DisplayName=sprintf('%s', formulation_names{formulation_idx}));
% 		% plot(0:num_nodes, one_sigma_minus, 'b-', DisplayName=sprintf('%s', formulation_names{formulation_idx}));
% 	end
% 	legend(legendUnq(), Location='northoutside', Orientation='horizontal', IconColumnWidth=15, FontSize=25)
% end
