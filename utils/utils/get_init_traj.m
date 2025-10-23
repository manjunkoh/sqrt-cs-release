function traj_nodes = get_init_traj(x_0, x_f, N, revs, MU, plot)
%GET_INIT_TRAJ Interpolates initial and final states using classical orbital elements
% FIXME: not great with handling revs 
arguments
    x_0 (6,1) double {mustBeNumeric}
    x_f (6,1) double {mustBeNumeric}
    N (1,1) double {mustBeInteger}
    revs (1,1) double {mustBeInteger} % revs: number of complete revolutions
    MU (1,1) double {mustBeNumeric}
    plot (1,1) logical = false
end

if x_0(3) == 0 && x_0 (6) == 0 & x_f(3) == 0 & x_f(6) == 0
    coe02D = sv2coe2D(x_0([1 2 4 5]), MU);
    coef2D = sv2coe2D(x_f([1 2 4 5]), MU);
    coef2D(end) = coef2D(end) + 2*pi*revs;
    coe2D = [
        linspace(coe02D(1), coef2D(1), N+1)
        linspace(coe02D(2), coef2D(2), N+1)
        linspace(coe02D(3), coef2D(3), N+1)
        linspace(coe02D(4), coef2D(4), N+1)
    ];
    traj_nodes = zeros(6,N+1);
    for i = 1:N+1
        traj_nodes([1 2 4 5],i) = coe2sv2D(coe2D(:,i), MU);
    end

else
    coe_0 = sv2coe(x_0, MU);
    coe_f = sv2coe(x_f, MU);
    % if coe_f(6) < coe_0(6)
    %     coe_f(6) = coe_f(6) + 2*pi;
    % end
    coe_f(6) = coe_f(6) + 2*pi*revs;
    coe = [
        linspace(coe_0(1), coe_f(1), N+1)
        linspace(coe_0(2), coe_f(2), N+1)
        linspace(coe_0(3), coe_f(3), N+1)
        linspace(coe_0(4), coe_f(4), N+1)
        linspace(coe_0(5), coe_f(5), N+1)
        linspace(coe_0(6), coe_f(6), N+1)
    ];
    traj_nodes = zeros(6,N+1); 
    for i = 1:N+1
        traj_nodes(:,i) = coe2sv(coe(:,i), MU);
    end
end

if plot
    figure
    plot3(traj_nodes(1,:), traj_nodes(2,:), traj_nodes(3,:), 'b')
    hold on
    plot3(x_0(1), x_0(2), x_0(3), 'r*')
    plot3(x_f(1), x_f(2), x_f(3), 'r*')
    grid on
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title('Initial Trajectory')
end
end