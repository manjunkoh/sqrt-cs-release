
function tl = plot_T_direct(tl, t_his, T, Tmax, varargin)

nexttile(tl,1)
hold on
stairs(t_his, [T(1, :) T(1,end)], varargin{:})
% xlim([0 t_his(end)])
ylabel('$T_x$')

nexttile(tl,2)
hold on
stairs(t_his, [T(2, :) T(2,end)], varargin{:})
ylabel('$T_y$')
% xlim([0 t_his(end)])

nexttile(tl,3)
hold on
stairs(t_his, [T(3, :) T(3,end)], varargin{:})
ylabel('$T_z$')
% xlim([0 t_his(end)])

nexttile(tl,4)
hold on
stairs(t_his, [norms(T) norm(T(:,end))], varargin{:})
yline(Tmax, 'k--', 'LineWidth', 1, 'HandleVisibility', 'Off')
text(t_his(end), Tmax, '$T_{\max}$', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top')
ylim([0 Tmax])
% xlim([0 t_his(end)])

xlabel('$t$ (nondimensional)')
ylabel('$\|T\|$')

end