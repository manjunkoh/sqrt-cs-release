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