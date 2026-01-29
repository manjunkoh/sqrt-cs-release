function add_zoomed_axis(ax, zoom_region, position, options)
    % Adds a secondary axis that zooms into a specified region of the original axis
    %
    % Parameters:
    % ax: Handle to the original axis
    % zoom_region: [x_min, x_max, y_min, y_max] specifying the region to zoom into
    % position: [x, y, width, height] specifying the position of the zoomed axis
    %   - If position_mode is 'figure' (default): position is in normalized figure coordinates [0-1]
    %   - If position_mode is 'data': x, y, width, and height are all in data coordinates of the original axes
    % zoom_factor: Factor by which to enlarge the zoomed region (optional)
    % position_mode: 'figure' (default) or 'data' - specifies coordinate system for position

	arguments
		ax (1,1) {mustBeA(ax, 'matlab.graphics.axis.Axes')} % Handle to the original axis
		zoom_region (1,4) {mustBeNumeric} % [x_min, x_max, y_min, y_max] specifying the region to zoom into
		position (1,4) {mustBeNumeric} % [x, y, width, height] specifying the position of the zoomed axis
		options.zoom_factor (1,1) {mustBeNumeric} = 1 % Factor by which to enlarge the zoomed region (optional)
		options.position_mode {mustBeMember(options.position_mode, {'figure', 'data'})} = 'data' % Coordinate system for position parameter
	end

	zoom_factor = options.zoom_factor;
	position_mode = options.position_mode;

    % Convert position from data coordinates to figure coordinates if needed
    if strcmp(position_mode, 'data')
        % Get the original axis position in normalized figure coordinates
        ax_pos = get(ax, 'Position');
        
        % Extract position values in data coordinates
        x_data = position(1);
        y_data = position(2);
        width_data = position(3);  % width in data coordinates (x-axis units)
        height_data = position(4);  % height in data coordinates (y-axis units)
        
        % Convert x, y from data coordinates to normalized figure coordinates
        x_norm = ax_pos(1) + ax_pos(3) * (x_data - ax.XLim(1)) / diff(ax.XLim);
        y_norm = ax_pos(2) + ax_pos(4) * (y_data - ax.YLim(1)) / diff(ax.YLim);
        
        % Convert width and height from data coordinates to normalized figure coordinates
        width_norm = width_data / diff(ax.XLim) * ax_pos(3);
        height_norm = height_data / diff(ax.YLim) * ax_pos(4);
        
        position = [x_norm, y_norm, width_norm, height_norm];
    end

    % Create a new axis for the zoomed region
    zoom_ax = axes('Position', position);

    % Copy the data from the original axis to the zoomed axis
    copyobj(allchild(ax), zoom_ax);

    % Set the limits of the zoomed axis to the specified region
    xlim(zoom_ax, zoom_region(1:2));
    ylim(zoom_ax, zoom_region(3:4));

    % Adjust the zoom factor if specified
    if zoom_factor ~= 1
        zoom_x_center = mean(zoom_region(1:2));
        zoom_y_center = mean(zoom_region(3:4));
        zoom_x_range = diff(zoom_region(1:2)) * zoom_factor;
        zoom_y_range = diff(zoom_region(3:4)) * zoom_factor;
        xlim(zoom_ax, [zoom_x_center - zoom_x_range / 2, zoom_x_center + zoom_x_range / 2]);
        ylim(zoom_ax, [zoom_y_center - zoom_y_range / 2, zoom_y_center + zoom_y_range / 2]);
    end

    % Add a box around the zoomed region in the original axis
    rectangle('Position', [zoom_region(1), zoom_region(3), diff(zoom_region(1:2)), diff(zoom_region(3:4))], ...
                     'EdgeColor', 'k', 'LineStyle', '-', 'Parent', ax);

    % Set the zoomed axis to clip anything outside its limits
    set(zoom_ax, 'Clipping', 'on');
    
    % Get grid properties from original axis if available
    if strcmp(ax.XGrid, 'on') || strcmp(ax.YGrid, 'on')
        grid_color_base = ax.GridColor;
        grid_alpha = ax.GridAlpha;
        grid_linestyle = ax.GridLineStyle;
        % Apply alpha to color (if color is RGB, convert to RGBA)
        if length(grid_color_base) == 3
            grid_color = [grid_color_base, grid_alpha];
        else
            grid_color = grid_color_base;
        end
    else
        grid_color = [0.15 0.15 0.15 1];  % Default grid color with alpha
        grid_linestyle = '-';
    end
    
    % Get axis limits for drawing grid lines
    x_lim = xlim(zoom_ax);
    y_lim = ylim(zoom_ax);
    
    % Generate reasonable tick positions for grid lines
    % Use a reasonable number of grid lines (adjust as needed)
    n_x_grid = 5;
    n_y_grid = 5;
    x_grid = linspace(x_lim(1), x_lim(2), n_x_grid);
    y_grid = linspace(y_lim(1), y_lim(2), n_y_grid);
    
    % Draw vertical grid lines
    hold(zoom_ax, 'on');
    for i = 1:length(x_grid)
        line(zoom_ax, [x_grid(i), x_grid(i)], y_lim, ...
            'Color', grid_color, 'LineStyle', grid_linestyle, ...
            'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
    
    % Draw horizontal grid lines
    for i = 1:length(y_grid)
        line(zoom_ax, x_lim, [y_grid(i), y_grid(i)], ...
            'Color', grid_color, 'LineStyle', grid_linestyle, ...
            'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
    hold(zoom_ax, 'off');
    
    % hide the xticks and yticks
    set(zoom_ax, 'XTick', [], 'YTick', []);

    % Connect the center of the top edge of the original rectangle to the center of the bottom edge of the zoomed axis
    % Convert data coordinates to normalized figure coordinates
    ax_pos = get(ax, 'Position');
    fig = ancestor(ax, 'figure');
    x_norm = @(x) ax_pos(1) + ax_pos(3) * (x - ax.XLim(1)) / diff(ax.XLim);
    y_norm = @(y) ax_pos(2) + ax_pos(4) * (y - ax.YLim(1)) / diff(ax.YLim);

    % Center of the top edge of the original rectangle (in data coordinates)
    top_center_x = mean(zoom_region(1:2));
    top_center_y = zoom_region(4);  % y_max
    
    % Center of the bottom edge of the zoomed axis (already in normalized figure coordinates)
    bottom_center_x = position(1) + position(3) / 2;
    bottom_center_y = position(2);
    
    % Draw line from top center of original rectangle to bottom center of zoomed axis
    % annotation(fig, 'line', ...
    %            [x_norm(top_center_x), bottom_center_x], ...
    %            [y_norm(top_center_y), bottom_center_y], ...
    %            'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);

    % Draw line from top left corner of original rectangle to top left corner of zoomed axis
    annotation(fig, 'line', ...
               [x_norm(zoom_region(1)), x_norm(position(1))], ...
               [y_norm(zoom_region(4)), y_norm(position(4))], ...
               'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
end