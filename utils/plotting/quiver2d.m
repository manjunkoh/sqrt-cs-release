function quiver2d(x, u, varargin)
    quiver(x(1,:), x(2,:), u(1,:), u(2,:), varargin{:});
end