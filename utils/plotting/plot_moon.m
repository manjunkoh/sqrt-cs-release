function plot_moon(mu, LU, varargin)
    % plot the moon in the CR3BP rotating frame
    R_MOON = 1737.1 / LU;
    x = 1 - mu;
    y = 0;
    z = 0;
    [X, Y, Z] = sphere(100);
    surf(x + R_MOON*X, y + R_MOON*Y, z + R_MOON*Z, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.5, varargin{:});
end
