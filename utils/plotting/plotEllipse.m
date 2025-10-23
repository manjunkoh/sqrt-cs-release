function varargout = plotEllipse(a,b,c,theta,varargin)
    t = linspace(0, 2 * pi);

    x = c(1) + a*cos(t)*cos(theta) - b*sin(t)*sin(theta);
    y = c(2) + a*cos(t)*sin(theta) + b*sin(t)*cos(theta);

    [varargout{1:nargout}] = plot(x,y, varargin{:});
end