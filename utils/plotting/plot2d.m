function varargout = plot2d(x, varargin)
    [varargout{1:nargout}] = plot(x(1,:), x(2,:), varargin{:});
    axis equal
end