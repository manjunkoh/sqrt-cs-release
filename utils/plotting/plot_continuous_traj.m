function plot_continuous_traj(z, SCvxProblem)

    [~, x, u, ~, ~, t] = getVarsFromz(z, SCvxProblem);

    f = SCvxProblem.DynamicalSystem.EoM;
    Nseg = SCvxProblem.TrajOpt.Nseg;
    x_arr = SCvxProblem.TrajOpt.x_arr;

    ode_options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    traj = ode45(@(t,x) f(t,x,u(:,1)), [t(1), t(2)], x(:,1), ode_options);

    for k = 2:Nseg
        traj = odextend(traj, @(t,x) f(t,x,u(:,k)), t(k+1), x(:,k), ode_options);
    end

    hold on
    plot(traj.y(1,:), traj.y(2,:))
    plot(x(1,:), x(2,:), 'o')
    plot(x_arr(1,:), x_arr(2,:), 'o')
    axis equal
end