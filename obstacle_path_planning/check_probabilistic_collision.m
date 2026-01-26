function [wall_safe_flag, obstacle_safe_flag, obstacle_safe_flags, constraint_violations] = check_probabilistic_collision(mu, P, obstacle_centers, obstacle_radii, wall_y_pos, state_risk)

    constraint_violation_tolerance = 1E-3; 

    if isempty(mu) || isempty(P)
        wall_safe_flag = NaN;
        obstacle_safe_flag = NaN;
        obstacle_safe_flags = [];
        constraint_violations = [];
        return;
    end
    wall_safe_flag = true;
    obstacle_safe_flags = true(size(obstacle_centers, 1), size(mu, 2));
    constraint_violations = zeros(size(obstacle_centers, 1), size(mu, 2));
    
    pos_idx = 1:2;
    z = norminv(1 - state_risk);
    for k = 1:size(mu, 2)
        for i = 1:size(obstacle_centers, 1)
            % consider a hyperplane that passes through the obstacle center and is perpendicular to the x-axis
            [a, b] = hyperplane_from_circular_obstacle(obstacle_centers(i,:)', obstacle_radii(i), mu(pos_idx,k));
            constraint_violations(i,k) = a' * mu(pos_idx,k) + b + z * sqrt(a' * P(pos_idx,pos_idx,k) * a);
            obstacle_safe_flags(i,k) = constraint_violations(i,k) <= constraint_violation_tolerance;
        end

        wall_safe_flag = wall_safe_flag && (mu(2,k) - wall_y_pos + z * sqrt(P(2,2,k)) <= 0);
    end

    constraint_violations = max(0, constraint_violations);

    obstacle_safe_flag = all(obstacle_safe_flags, 'all');

    return

    figure(Position=[0, 0, 15, 15]);
    hold on

    plot_obstacles(obstacle_centers, obstacle_radii);
    plot_wall(wall_y_pos);

    axis equal

    % xlim([-1.4, 11])
    % ylim([-3, 2.0])

    for k = 1:size(mu, 2)
        for i = 1:size(obstacle_centers, 1)
            % if constraint_violations(i,k) > 0
                if constraint_violations(i,k) > 0
                    FaceColor = 'r';
                else
                    FaceColor = '#0082B2';
                end

                fill_confidence_ellipse(mu(pos_idx,k), P(pos_idx,pos_idx,k), 1-state_risk, '', FaceColor=FaceColor, FaceAlpha=0.5, EdgeColor='none', DisplayName="$3 \sigma$ ellipse");

                [a, b] = hyperplane_from_circular_obstacle(obstacle_centers(i,:)', obstacle_radii(i), mu(pos_idx,k));
                constraint_violation = a' * mu(pos_idx,k) + b + z * sqrt(a' * P(pos_idx, pos_idx, k) * a);
                % plot_hyperplane(a, b, [3, 7], 'k-', LineWidth=1);
                % line([obstacle_centers(i,1), mu(1,k)], [obstacle_centers(i,2), mu(2,k)], 'Color', 'r', 'LineWidth', 1);
                drawnow
                keyboard
            % end

        end
    end

end
