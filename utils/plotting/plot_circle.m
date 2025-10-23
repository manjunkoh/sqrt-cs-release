function varargout = plot_circle(center,rs, varargin)
    t = linspace(0,2*pi);
    hold on;
    for r = rs
        [varargout{1:nargout}] = plot(center(1)+r*cos(t), center(2) + r*sin(t), varargin{:});
    end
end
