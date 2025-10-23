function tl = plot_ZOH_control(t_his, u_his, tl, varargin)
    if nargin == 2
        figure
        tl = tiledlayout(4,1);
        tl.TileSpacing = 'tight';
        tl.Padding = 'tight';
    end

    u_his_stairs = [u_his u_his(:,end)];

    nexttile(tl, 1)
    stairs(t_his, u_his_stairs(1,:), varargin{:})
    hold on

    nexttile(tl, 2)
    stairs(t_his, u_his_stairs(2,:), varargin{:})
    hold on

    nexttile(tl, 3)
    stairs(t_his, u_his_stairs(3,:), varargin{:})
    hold on

    nexttile(tl, 4)
    stairs(t_his, norms(u_his_stairs(1:3,:)), varargin{:})
    hold on

end