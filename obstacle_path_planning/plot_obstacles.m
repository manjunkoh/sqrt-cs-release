function plot_obstacles(centers, radii)
    hold on
    for i = 1:size(centers, 1)
        filled_circle(centers(i, :), radii(i));
        text(centers(i, 1), centers(i, 2), sprintf('%d', i));
    end
end
