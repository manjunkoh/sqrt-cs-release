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
