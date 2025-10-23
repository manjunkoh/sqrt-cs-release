function plot_Earth(mu, LU)
    % plot the Earth in the CR3BP rotating frame
    R_EARTH = 6371 / LU;
    x = -mu;
    y = 0;
    z = 0;
    [X, Y, Z] = sphere(100);
    surf(x + R_EARTH*X, y + R_EARTH*Y, z + R_EARTH*Z, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end