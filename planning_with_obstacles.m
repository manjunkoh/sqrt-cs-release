clc; clear;

%%
function plot_obstacles(centers, radius)
    hold on
    for i = 1:size(centers, 1)
        filled_circle(centers(i, :), radius);
    end
end

function filled_circle(center, radius)
    rectangle('Position', [center(1)-radius, center(2)-radius, 2*radius, 2*radius], 'Curvature', [1, 1], 'FaceColor', 'k', 'FaceAlpha', 0.5);
end

function plot_wall(wall_y_pos)
    yline(wall_y_pos, 'k--');
end

function plot_problem(obstacle_center, obstacle_radius, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, options)
    arguments
        obstacle_center
        obstacle_radius
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
    plot_obstacles(obstacle_center, obstacle_radius);
    % plot_wall(wall_y_pos);
    plot3sigmaEllipse(mu_0, Sigma_0, Color='#D55E00', DisplayName='Start')
    plot3sigmaEllipse(mu_f, Sigma_f, Color='#D55E00', LineStyle=":", DisplayName='Goal')
    xlabel('$x$')
    ylabel('$y$')
end

wall_y_pos = 2.2;
obstacle_center = [5, 0];
obstacle_radius = 1.2;

% plot_problem(obstacle_center, obstacle_radius, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f);

num_nodes = 50;
t_f = 30;
dt = t_f / num_nodes;
nx = 4;
nu = 2;

mu_0 = [0; 0; 0; 0];
mu_f = [10; 0; 0; 0];

Sigma_0 = diag([0.1, 0.1, 0.01, 0.01]);
Sigma_f = Sigma_0;

Q = 0.001 * eye(nx);
R = 0.01 * eye(nu);

A = [eye(2), dt*eye(2); zeros(2, 2), eye(2)];
B = [0.5*dt^2*eye(2); dt*eye(2)];
% G = 0.01 * [zeros(2, 2); eye(2)];
G = [0.01 * eye(2); 0.01 * eye(2)];
% q = 0.1;
% G = sqrt(q * dt) * eye(4);

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

circular_obstacles={struct('center', obstacle_center, 'radius', obstacle_radius, 'p', state_risk)};

%% Solve the obstacle-free deterministic problem
yalmip('clear')
x = sdpvar(nx, num_nodes+1, 'full');
u = sdpvar(nu, num_nodes, 'full');

constraints = [];
objective = 0;

for k = 1:num_nodes
    constraints = [constraints, x(:,k+1) == A * x(:,k) + B * u(:,k)];
    % constraints = [constraints, x(2,k) <= wall_y_pos];
    constraints = [constraints, u(:,k) <= u_max];
    constraints = [constraints, u(:,k) >= -u_max];
end

constraints = [constraints, x(:,1) == mu_0, x(:,num_nodes+1) == mu_f];

for k = 1:num_nodes
    objective = objective + x(:,k)' * Q * x(:,k) + u(:,k)' * R * u(:,k);
end

sol = optimize(constraints, objective);

if sol.problem
    disp('Unconstrained problem infeasible');
    return
end

x_opt = value(x);
u_opt = value(u);

plot_problem(obstacle_center, obstacle_radius, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f);
plot(x_opt(1,:), x_opt(2,:), 'r.-');

%% Solve the constrained deterministic problem via iterative approach
% Linearize the obstacles around the current state
penalty_scalar = 10;
max_iters = 50;

yalmip('clear')
x = sdpvar(nx, num_nodes+1, 'full');
u = sdpvar(nu, num_nodes, 'full');
lambda = sdpvar(1, num_nodes, 'full');

x_ref = x_opt;

for iter = 1:max_iters

    fprintf("Iter %d\n", iter)

    constraints = [];
    objective = 0;

    for k = 1:num_nodes
        constraints = [constraints, x(:,k+1) == A * x(:,k) + B * u(:,k)];
        constraints = [constraints, u(:,k) <= u_max];
        constraints = [constraints, u(:,k) >= -u_max];
    end
    if iter == 1
        constraints = [constraints, x(2, ceil(num_nodes* 0.5)) >= 0.1]; % symmetry-breaking constraint
    end
    constraints = [constraints, x(:,1) == mu_0, x(:,num_nodes+1) == mu_f];

    for k = 1:num_nodes
        a = - (x_ref(1:2,k) - obstacle_center');
        b = - 0.5 * norm(x_ref(1:2,k) - obstacle_center')^2  + 0.5 * obstacle_radius^2 - a' * x_ref(1:2,k);
        constraints = [constraints
            a' * x(1:2,k) + b <= lambda(k)
        ];
    end
    
    constraints = [constraints, lambda >= 0];

    for k = 1:num_nodes
        objective = objective + x(:,k)' * Q * x(:,k) + u(:,k)' * R * u(:,k) + penalty_scalar * lambda(k);
    end

    sol = optimize(constraints, objective, sdpsettings('verbose', 0));

    if sol.problem
        disp('Constrained problem infeasible');
        break;
    end

    x_opt = value(x);
    u_opt = value(u);

    if norm(x_opt - x_ref) < 1e-3 && is_collision_free(x_opt, obstacle_center, obstacle_radius)

        fprintf('Converged in %d iterations\n', iter);
        break;
    end

    x_ref = x_opt;

    if iter == max_iters
        fprintf('Reached max iters.\n')
        return
    end

end

function yn = is_collision_free(x_opt, obstacle_center, obstacle_radius)
    yn = false;
    for k = 1:size(x_opt, 2)
        if norm(x_opt(1:2,k) - obstacle_center') < obstacle_radius
            return
        end
    end
    yn = true;
end
        
%% Plot the results
plot_problem(obstacle_center, obstacle_radius, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f);
plot(x_opt(1,:), x_opt(2,:), 'r.-');

%% Solve the constrained stochastic problem via iterative approach
max_iters = 100;
% penalty_scalar = 100;
convergence_tolerance = 1E-2;

yalmip('clear')
mu = sdpvar(nx, num_nodes+1, 'full');
v = sdpvar(nu, num_nodes, 'full');
P = sdpvar(nx, nx, num_nodes+1);
U = sdpvar(nu, nx, num_nodes, 'full');
Y = sdpvar(nu, nu, num_nodes);
% lambda = sdpvar(1, num_nodes, 'full');

mu_ref = x_opt;
P_ref = repmat(Sigma_0, [1,1,num_nodes+1]);
v_ref = u_opt;
Y_ref = repmat(0.00001 * eye(nu), [1,1,num_nodes]);

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

    % for k = 1:num_nodes
    %     sqrt_ref = sqrt(P_ref(2,2,k));
    %     constraints = [constraints,
    %         z / (2 * sqrt_ref) * (P(2,2,k)) + z * sqrt_ref / 2 + mu(2,k) - wall_y_pos <= 0
    %     ];
    % end

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
                z / (2 * sqrt_ref) * (a' * Y_k * a) + a' * v_k - b + z * sqrt_ref / 2 <= 0
            ];
        end
    end

    % State chance constraints
    for k = 1:num_nodes+1
        a = - (x_ref(1:2,k) - obstacle_center');
        b = - 0.5 * norm(x_ref(1:2,k) - obstacle_center')^2 + 0.5 * obstacle_radius^2 - a' * x_ref(1:2,k);
        P_ref_pos_k = P_ref(1:2,1:2,k);
        sqrt_ref = sqrt(a' * P_ref_pos_k * a);
        if sqrt_ref <= 0
            % numerically, the covariance can be non-PSD; in this case
            % simply set this value to a small positive number
            % Using nearestSPD sometimes does not terminate for a long
            % time, so it isn't used here
            sqrt_ref = 0.001;
        end
        constraints = [constraints
            z / (2 * sqrt_ref) * (a' * P(1:2,1:2,k) * a) + z * sqrt_ref / 2 + a' * mu(1:2,k) + b <= 0
        ];
    end

    % constraints = [constraints, lambda >= 0];

    for k = 1:num_nodes
        objective = objective + mu(:,k)' * Q * mu(:,k) + v(:,k)' * R * v(:,k) + trace(Q * P(:,:,k)) + trace(R * Y(:,:,k));
    end

    % objective_augmented = objective + penalty_scalar * sum(lambda);

    fprintf("Iteration %d    ", iter)

    sol = optimize(constraints, objective, sdpsettings('verbose', 0));


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
    % lambda_opt = value(lambda);
    
    objective_full_covariance = value(objective);

    if all(vecnorm(mu_opt - mu_ref, Inf) < convergence_tolerance)
        time_full_covariance = toc;
        fprintf('Converged in %d iterations in %.3f seconds\n', iter, time_full_covariance);
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

%% Plot the results
plot_problem(obstacle_center, obstacle_radius, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f);
for k = 1:num_nodes
    fill3sigmaEllipse(mu_opt(:,k), P_opt(1:2,1:2,k), '', FaceColor='#0082B2', FaceAlpha=0.5, EdgeColor='none', DisplayName="$3 \sigma$ ellipse");
end
plot(mu_opt(1,:), mu_opt(2,:), 'k.-', DisplayName='mean', LineWidth=1);
legend(legendUnq(), Location='northoutside', Orientation='horizontal', IconColumnWidth=15, FontSize=25)

xlim([-1.4, 11])
ylim([-3, 1.5])
exportgraphics(gcf, 'figures/planning_with_obstacles_full_covariance.png', Resolution=300)
exportgraphics(gcf, 'figures/planning_with_obstacles_full_covariance.pdf', ContentType='vector')

%% Solve with SQRT QR method
yalmip('clear')
init_guess = struct();
init_guess.S = interpolate_lower_triangular(chol(Sigma_0, 'lower'), chol(Sigma_f, 'lower'), num_nodes+1, 'log-cholesky');
init_guess.L = zeros(nu, nx, num_nodes);
init_guess.mu = x_opt;
init_guess.v = u_opt;

% chance_constraints_state={struct('type', 'affine', 'alpha', ...
    % [0; 1; 0; 0], 'beta', wall_y_pos, 'p', state_risk, 'nodes', 1:num_nodes+1)};

prob_qr = SqrtQRCovarianceSteering(init_guess, ...
    N=num_nodes, ...
    A_sys=A_sys, B_sys=B_sys, G_sys=G_sys, ...
    objective_type='LQR', ...
    mu_0=mu_0, mu_f=mu_f, ...
    P_0=Sigma_0, P_f=Sigma_f, Q=Q, R=R, ...
    chance_constraints_control=chance_constraints_control,...
    circular_obstacles=circular_obstacles);

scp_params = SCPParams();
scp_params.k_max = 100;
scp_params.tol_opt = 1E-2;
scp_params.tol_feas = 1E-4;
scp_params.w_init = 100;
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

%% Plot the results
plot_problem(obstacle_center, obstacle_radius, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f);
for k = 1:num_nodes
    fill3sigmaEllipse(prob_qr.mu(:,k), prob_qr.P(:,:,k), '', FaceColor='#0082B2', FaceAlpha=0.5, EdgeColor='none', DisplayName="$3 \sigma$ ellipse");
end
plot(prob_qr.mu(1,:), prob_qr.mu(2,:), 'k.-', DisplayName='mean', LineWidth=1);

legend(legendUnq(), Location='northoutside', Orientation='horizontal', IconColumnWidth=15, FontSize=25)

xlim([-1.4, 11])
ylim([-3, 1.5])
exportgraphics(gcf, 'figures/planning_with_obstacles_sqrt_qr.png', Resolution=300)
exportgraphics(gcf, 'figures/planning_with_obstacles_sqrt_qr.pdf', ContentType='vector')

%% Plot the results for both methods using a single figure
% figure;
% tiledlayout(2, 1);
% nexttile;
% 
% plot_problem(obstacle_center, obstacle_radius, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);
% for k = 1:num_nodes
%     fill3sigmaEllipse(mu_opt(:,k), P_opt(1:2,1:2,k), '', FaceColor='#0082B2', FaceAlpha=0.5, EdgeColor='none', DisplayName="$3 \sigma$ ellipse");
% end
% plot(mu_opt(1,:), mu_opt(2,:), 'k.-', DisplayName='mean', LineWidth=1);
% xlim([-2, 11])
% ylim([-3, 1.5])
% nexttile;
% 
% plot_problem(obstacle_center, obstacle_radius, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, fig=gcf);
% for k = 1:num_nodes
%     fill3sigmaEllipse(prob_qr.mu(:,k), prob_qr.P(:,:,k), '', FaceColor='#0082B2', FaceAlpha=0.5, EdgeColor='none', DisplayName="$3 \sigma$ ellipse");
% end
% plot(prob_qr.mu(1,:), prob_qr.mu(2,:), 'k.-', DisplayName='mean', LineWidth=1);
% 
% legend(legendUnq(), Location='south', Orientation='horizontal', IconColumnWidth=15)
% xlim([-2, 11])
% ylim([-3, 1.5])
%% Display summary table
fprintf('\n');
fprintf('================================================================================\n');
fprintf('SUMMARY TABLE\n');
fprintf('================================================================================\n');
fprintf('%-30s | %-12s | %-15s | %-15s\n', 'Method', 'Iterations', 'Cost', 'Time (s)');
fprintf('--------------------------------------------------------------------------------\n');

% Format full covariance row
if isnan(iters_full_covariance) || isnan(objective_full_covariance) || isnan(time_full_covariance)
    fprintf('%-30s | %-12s | %-15s | %-15s\n', 'Full Covariance', 'N/A', 'N/A', 'N/A');
else
    fprintf('%-30s | %-12d | %-15.6f | %-15.3f\n', 'Full Covariance', iters_full_covariance, objective_full_covariance, time_full_covariance);
end

% Format sqrt QR row
if isnan(iters_qr) || isnan(objective_qr) || isnan(time_qr)
    fprintf('%-30s | %-12s | %-15s | %-15s\n', 'Sqrt QR', 'N/A', 'N/A', 'N/A');
else
    fprintf('%-30s | %-12d | %-15.6f | %-15.3f\n', 'Sqrt QR', iters_qr, objective_qr, time_qr);
end

fprintf('================================================================================\n');
fprintf('\n');

%% Additional plotting for the sqrt method
%% Plot the coordinate-wise control 
figure;
tiledlayout(nu, 1);
for j = 1:nu
    nexttile;
    plot(prob_qr.v(j,:), 'b.-');
    hold on;
    two_sigma_upper = prob_qr.v(j,:) + 2 * sqrt(squeeze(prob_qr.P_u(j,j,:))');
    two_sigma_lower = prob_qr.v(j,:) - 2 * sqrt(squeeze(prob_qr.P_u(j,j,:))');
    plot(two_sigma_upper, 'r--');
    plot(two_sigma_lower, 'r--');
    yline(u_max, 'k--')
    yline(-u_max, 'k--')
    ylabel(sprintf('$u_%d$', j));
end
xlabel('Node');

%% Plot the coordinate-wise state
figure;
tiledlayout(nx, 1);
for i = 1:nx
    nexttile;
    plot(prob_qr.mu(i,:), 'b.-');
    hold on;
    two_sigma_upper = prob_qr.mu(i,:) + 2 * sqrt(squeeze(prob_qr.P(i,i,:))');
    two_sigma_lower = prob_qr.mu(i,:) - 2 * sqrt(squeeze(prob_qr.P(i,i,:))');
    plot(two_sigma_upper, 'r--');
    plot(two_sigma_lower, 'r--');
    ylabel(sprintf('$x_%d$', i));
end
xlabel('Node');