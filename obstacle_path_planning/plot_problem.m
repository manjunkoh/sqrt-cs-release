function plot_problem(obstacle_center, obstacle_radii, wall_y_pos, mu_0, Sigma_0, mu_f, Sigma_f, options)
    arguments
        obstacle_center
        obstacle_radii
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
    plot_obstacles(obstacle_center, obstacle_radii);
    plot_wall(wall_y_pos);
    plot3sigmaEllipse(mu_0, Sigma_0, Color='#D55E00', DisplayName='Start')
    plot3sigmaEllipse(mu_f, Sigma_f, Color='#D55E00', LineStyle=":", DisplayName='Goal')
    xlabel('$x$')
    ylabel('$y$')
end
