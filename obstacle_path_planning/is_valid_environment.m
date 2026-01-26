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
