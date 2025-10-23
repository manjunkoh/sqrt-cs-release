function varargout = color_along_line(varargin, options)

    arguments (Repeating)
        varargin
    end

    arguments
        options.LineWidth = 3
    end

    if nargin < 3
        error('Not enough input arguments.')
    elseif nargin > 4
        error('Too many input arguments.')
    elseif nargin == 3
        x = varargin{1}; x = x(:)';
        y = varargin{2}; y = y(:)';
        data = varargin{3}; data = data(:)';
        p = patch([x nan] , [y nan], [data nan], 'EdgeColor', 'interp', 'FaceColor', 'none', 'LineWidth', options.LineWidth, HandleVisibility='off');
    else
        x = varargin{1}; x = x(:)';
        y = varargin{2}; y = y(:)';
        z = varargin{3}; z = z(:)';
        data = varargin{4}; data = data(:)';
        p = patch([x nan] , [y nan], [z nan], [data nan], 'EdgeColor', 'interp', 'FaceColor', 'none', 'LineWidth', 3);
    end

    if nargout > 0
        [varargout{1}] = colorbar();
        varargout{2} = p;
    end
end    