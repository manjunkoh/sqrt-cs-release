function add_zoomed_axis(ax, zoom_region, position, options)
    % Adds a secondary axis that zooms into a specified region of the original axis
    %
    % Parameters:
    % ax: Handle to the original axis
    % zoom_region: [x_min, x_max, y_min, y_max] specifying the region to zoom into
    % position: [x, y, width, height] specifying the position of the zoomed axis
    % zoom_factor: Factor by which to enlarge the zoomed region (optional)
    % vertices_to_connect: Indices of rectangle vertices to connect to the zoomed axis (optional)

	arguments
		ax (1,1) {mustBeA(ax, 'matlab.graphics.axis.Axes')} % Handle to the original axis
		zoom_region (1,4) {mustBeNumeric} % [x_min, x_max, y_min, y_max] specifying the region to zoom into
		position (1,4) {mustBeNumeric} % [x, y, width, height] specifying the position of the zoomed axis
		options.zoom_factor (1,1) {mustBeNumeric} = 1 % Factor by which to enlarge the zoomed region (optional)
		options.vertices_to_connect (1,:) {mustBeNumeric} = 1:4 % Indices of rectangle vertices to connect to the zoomed axis (optional, [bottom-left, bottom-right, top-left, top-right])
	end

	zoom_factor = options.zoom_factor;
	vertices_to_connect = options.vertices_to_connect;

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
                     'EdgeColor', 'r', 'LineStyle', '--', 'Parent', ax);

    % Set the zoomed axis to clip anything outside its limits
    set(zoom_ax, 'Clipping', 'on');

    % Connect the specified edges of the rectangle to the zoomed axis
    % Convert data coordinates to normalized figure coordinates
    ax_pos = get(ax, 'Position');
    fig = ancestor(ax, 'figure');
    x_norm = @(x) ax_pos(1) + ax_pos(3) * (x - ax.XLim(1)) / diff(ax.XLim);
    y_norm = @(y) ax_pos(2) + ax_pos(4) * (y - ax.YLim(1)) / diff(ax.YLim);

    rect_vertices = [zoom_region(1), zoom_region(3); ... % Bottom-left
                     zoom_region(1) + diff(zoom_region(1:2)), zoom_region(3); ... % Bottom-right
                     zoom_region(1), zoom_region(3) + diff(zoom_region(3:4)); ... % Top-left
                     zoom_region(1) + diff(zoom_region(1:2)), zoom_region(3) + diff(zoom_region(3:4))]; % Top-right

    for i = vertices_to_connect
        annotation(fig, 'line', ...
                   [x_norm(rect_vertices(i, 1)), position(1) + position(3) * (i == 2 || i == 4)], ...
                   [y_norm(rect_vertices(i, 2)), position(2) + position(4) * (i > 2)], ...
                   'Color', 'k', 'LineStyle', '-', 'LineWidth', 1);
    end
end