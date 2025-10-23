function plot_traj_w_thrust(traj, u_handle, throttle_threshold)

    arguments
        traj struct
        u_handle function_handle
        throttle_threshold = 0.1
    end

    % Extract the state and control history
    t = linspace(traj.x(1), traj.x(end), 10000);
    x = deval(traj, t);

    nu = length(u_handle(traj.x(1)));
    u = zeros(nu, length(t));
    for i = 1:length(t)
        u(:,i) = u_handle(t(i));
    end

    u_norms = norms(u);
    u_max = max(u_norms);

    % plot the trajectory in red if the thrust is larger than a specified 
    % percentage of the max thrust

    traj_thrust = x;
    traj_coast = x;

    traj_thrust(:, u_norms < throttle_threshold*u_max) = NaN;
    traj_coast(:, u_norms >= throttle_threshold*u_max) = NaN;

    hold on;
    plot3d(traj_thrust, 'r');
    plot3d(traj_coast, 'k');

    % plot the trajectory along with the thrust magnitude over time
    % u_norms = norms(u);
    % u_norms = [u_norms u_norms(end)];
    
    % surface([x(1,:); x(1,:)],...
    %         [x(2,:); x(2,:)],...
    %         [x(3,:); x(3,:)],...
    %         [u_norms;u_norms],...
    %         'facecol','no',...
    %         'edgecol','interp',...
    %         'linew',2);
    % axis equal
    % colormap('autumn')
    % colorbar
end