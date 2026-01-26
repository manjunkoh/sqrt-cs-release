function [collision_flags_all] = count_collision_samples(x_hist_all, obstacle_centers, obstacle_radii)
    collision_flags_all = false(size(x_hist_all, 2), size(obstacle_centers, 1), size(x_hist_all, 3));
    for i = 1:size(x_hist_all, 3)
        x_hist = x_hist_all(:,:,i);
        collision_flags = ~is_collision_free_per_node_per_obstacle(x_hist, obstacle_centers, obstacle_radii);
        collision_flags_all(:,:,i) = collision_flags;
    end

    collision_flags_all = sum(collision_flags_all, 3);
end
