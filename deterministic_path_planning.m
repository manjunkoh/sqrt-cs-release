%% Solve the obstacle-free deterministic problem
% yalmip('clear')
% x = sdpvar(nx, num_nodes+1, 'full');
% u = sdpvar(nu, num_nodes, 'full');
% 
% constraints = [];
% objective = 0;
% 
% for k = 1:num_nodes
%     constraints = [constraints, x(:,k+1) == A * x(:,k) + B * u(:,k)];
%     % constraints = [constraints, x(2,k) <= wall_y_pos];
%     constraints = [constraints, u(:,k) <= u_max];
%     constraints = [constraints, u(:,k) >= -u_max];
% end
% 
% constraints = [constraints, x(:,1) == mu_0, x(:,num_nodes+1) == mu_f];
% 
% for k = 1:num_nodes
%     objective = objective + x(:,k)' * Q * x(:,k) + u(:,k)' * R * u(:,k);
% end
% 
% sol = optimize(constraints, objective);
% 
% if sol.problem
%     disp('Unconstrained problem infeasible');
%     return
% end
% 
% x_opt = value(x);
% u_opt = value(u);
% 
% plot_problem(obstacle_centers, obstacle_radius, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f);
% plot(x_opt(1,:), x_opt(2,:), 'r.-');

%% Solve the constrained deterministic problem via iterative approach
% Linearize the obstacles around the current state
% penalty_scalar_u = 10;
% max_iters = 50;
% 
% yalmip('clear')
% x = sdpvar(nx, num_nodes+1, 'full');
% u = sdpvar(nu, num_nodes, 'full');
% lambda = sdpvar(1, num_nodes, 'full');
% 
% % x_ref = x_opt;
% x_ref = linspace_vec(mu_0, mu_f, num_nodes+1);
% 
% for iter = 1:max_iters
% 
%     fprintf("Iter %d\n", iter)
% 
%     constraints = [];
%     objective = 0;
% 
%     for k = 1:num_nodes
%         constraints = [constraints, x(:,k+1) == A * x(:,k) + B * u(:,k)];
%         constraints = [constraints, u(:,k) <= u_max];
%         constraints = [constraints, u(:,k) >= -u_max];
%     end
%     % if iter == 1
%         % constraints = [constraints, x(2, ceil(num_nodes* 0.5)) >= 0.1]; % symmetry-breaking constraint
%     % end
%     constraints = [constraints, x(:,1) == mu_0, x(:,num_nodes+1) == mu_f];
% 
%     for k = 1:num_nodes
%         for i = 1:size(obstacle_centers, 1)
%             a = - (x_ref(1:2,k) - obstacle_centers(i,:)');
%             b = - 0.5 * norm(x_ref(1:2,k) - obstacle_centers(i,:)')^2  + 0.5 * obstacle_radius^2 - a' * x_ref(1:2,k);
%             constraints = [constraints
%                 a' * x(1:2,k) + b <= lambda(k)
%             ];
%         end
%     end
% 
%     constraints = [constraints, lambda >= 0];
% 
%     for k = 1:num_nodes
%         objective = objective + x(:,k)' * Q * x(:,k) + u(:,k)' * R * u(:,k) + penalty_scalar_u * lambda(k);
%     end
% 
%     sol = optimize(constraints, objective, sdpsettings('verbose', 0));
% 
%     if sol.problem
%         disp('Constrained problem infeasible');
%         break;
%     end
% 
%     x_opt = value(x);
%     u_opt = value(u);
% 
%     if norm(x_opt - x_ref) < 1e-3 && is_collision_free(x_opt, obstacle_centers, obstacle_radius)
% 
%         fprintf('Converged in %d iterations\n', iter);
%         break;
%     end
% 
%     x_ref = x_opt;
% 
%     if iter == max_iters
%         fprintf('Reached max iters.\n')
%         return
%     end
% 
% end
  
%% Plot the results
% plot_problem(obstacle_centers, obstacle_radius, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f);
% plot(x_opt(1,:), x_opt(2,:), 'r.-');