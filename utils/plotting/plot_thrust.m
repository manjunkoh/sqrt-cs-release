function plot_thrust(z, Problem)

    N = Problem.TrajOpt.Nseg;
    nx = Problem.TrajOpt.Nst;
    nu = Problem.TrajOpt.Nct;
    [~, x, u] = getVarsFromz(z, Problem);
    t_his = Problem.TrajOpt.t_his;
    log_mass = x(end,:);
    mass = exp(log_mass);
    T = u(1:nu - 1, :) .* mass(1:end - 1);

    tl = tiledlayout(4, 1);
    tl.TileSpacing = 'tight';
    tl.Padding = 'tight';

    nexttile
    stairs(t_his(1:end), [T(1, :) T(1,end)], 'k-')
    ylabel('$T_x$')

    nexttile
    stairs(t_his(1:end), [T(2, :) T(2,end)], 'k-')
    ylabel('$T_y$')

    nexttile
    stairs(t_his(1:end), [T(3, :) T(3,end)], 'k-')
    ylabel('$T_z$')

    nexttile
    hold on
    stairs(t_his(1:end), [norms(T) norm(T(:,end))], 'k-')
    yline(Problem.TrajOpt.Tmax, 'r--', 'LineWidth', 2)

    xlabel('$t$ (non-dim)')
    ylabel('$\|T\|$')