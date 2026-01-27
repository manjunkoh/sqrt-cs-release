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