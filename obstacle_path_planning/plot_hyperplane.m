function plot_hyperplane(a, b, xlim, varargin)
    x = linspace(xlim(1), xlim(2), 100);
    y = (-a(1) * x - b) / a(2);
    plot(x, y, varargin{:});
end
