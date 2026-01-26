function flags = is_collision_free_per_node_per_obstacle(x_hist, obstacle_centers, obstacle_radii)
    flags = true(size(x_hist, 2), size(obstacle_centers, 1));
    for k = 1:size(x_hist, 2)
        for i = 1:size(obstacle_centers, 1)
            if norm(x_hist(1:2,k) - obstacle_centers(i,:)') < obstacle_radii(i)
                flags(k,i) = false;
            end
        end
    end
end
