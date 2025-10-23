
function tl = plot_T_indirect(tl, t, T, Tmax, varargin)

nexttile(tl,1)
hold on
plot(t, T(1,:), varargin{:})
ylabel('$T_x$')

nexttile(tl,2)
hold on
plot(t, T(2,:), varargin{:})
ylabel('$T_y$')

nexttile(tl,3)
hold on
plot(t, T(3,:), varargin{:})
ylabel('$T_z$')

nexttile(tl,4)
hold on
plot(t, norms(T), varargin{:})
yline(Tmax, 'k--', 'LineWidth', 1, 'handleVisibility', 'Off')

xlabel('$t$ (nondimensional)')
ylabel('$\|T\|$')

end