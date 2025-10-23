function plot_control_from_z(z, Problem)

    N = Problem.TrajOpt.Nseg;
    nx = Problem.TrajOpt.Nst;
    nu = Problem.TrajOpt.Nct;
    [x, u_his] = getVarsFromz(z, nx, nu, N);
    t_his = Problem.TrajOpt.t_his;

    tl = tiledlayout(4,1);
    tl.TileSpacing = 'tight';
    tl.Padding = 'tight';

    nexttile
    stairs(t_his(1:end-1), u_his(1,:), 'k-')

    nexttile
    stairs(t_his(1:end-1), u_his(2,:), 'k-')

    nexttile
    stairs(t_his(1:end-1), u_his(3,:), 'k-')

    nexttile
    stairs(t_his(1:end-1), norms(u_his(1:3,:)), 'k-')

