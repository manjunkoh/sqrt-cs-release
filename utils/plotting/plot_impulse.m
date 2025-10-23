function plot_impulse(t_his, u_his)
    % Plot the magnitude of the impulsive control input
    stem(t_his, vecnorm(u_his), 'filled', 'LineWidth', 1.5);
end