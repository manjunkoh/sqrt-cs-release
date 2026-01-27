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