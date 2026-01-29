function [a, b] = hyperplane_from_circular_obstacle(obstacle_center, obstacle_radius, pos)
    % Returns a hyperplane that is tangent to the circlular obstacle and is
    % perpendicular to the line between the obsacle center and pos
	c = obstacle_center;
	r = obstacle_radius;
	n = pos - c;
	a = - n;
	b = n' * c + norm(n) * r;
end