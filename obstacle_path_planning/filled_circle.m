function filled_circle(center, radius)
    rectangle('Position', [center(1)-radius, center(2)-radius, 2*radius, 2*radius], 'Curvature', [1, 1], 'FaceColor', 'k', 'FaceAlpha', 0.5);
end
