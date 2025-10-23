function traj = getcontinuoustraj(f, N, x0, t_his, u, odeoptions, ode)
    %GETCONTINUOUSTRAJ Get a continuous trajectory from the control solution
    arguments
        f (1,1) function_handle
        N (1,1) {mustBeInteger, mustBePositive}
        x0 (:,1) double
        t_his  double
        u double
        odeoptions struct = odeset('RelTol', 1e-12, 'AbsTol', 1e-12)
        ode = @ode45
    end
    % Get the continuous trajectory
    traj = ode(@(t,x) f(t,x,u(:,1)), [t_his(1) t_his(2)], x0, odeoptions);
    for k = 2:N
        traj = odextend(traj, @(t,x) f(t,x,u(:,k)),  t_his(k+1));
    end
end
