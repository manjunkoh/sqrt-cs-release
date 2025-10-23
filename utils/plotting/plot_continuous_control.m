function tl = plot_continuous_control(t, u, tl)
    if nargin == 2
        figure
        tl = tiledlayout(4,1);
        tl.TileSpacing = 'tight';
        tl.Padding = 'tight';
    end

    nexttile(tl, 1)
    plot(t, u(1,:), 'k-')
    hold on

    nexttile(tl, 2)
    plot(t, u(2,:), 'k-')
    hold on

    nexttile(tl, 3)
    plot(t, u(3,:), 'k-')
    hold on

    nexttile(tl, 4)
    plot(t, norms(u), 'k-')
    hold on