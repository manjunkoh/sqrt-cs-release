function figure_settings()

% text interpreter
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultConstantLineInterpreter', 'latex');
set(groot, 'defaultTextboxshapeInterpreter', 'latex');
set(groot, 'defaultTextarrowshapeInterpreter', 'latex');

set(groot, 'defaultTextboxshapeFaceAlpha', 0.5)

% GUI font
% set(groot, 'defaultUicontrolFontName', 'Times');

% Axis font
FontName = 'Arial';
set(groot, 'defaultAxesFontName', FontName);
set(groot, 'defaultTextFontName', FontName);
set(groot, 'defaultLegendFontName', FontName);
set(groot, 'defaultConstantLineFontName', FontName);

% GUI fontsize
set(groot, 'defaultUicontrolFontSize', 9);

% Set figure background color; for saving with export_fig
set(groot, 'defaultFigureColor', 'white')
% set(groot, 'defaultFigureColor', 'black')
set(groot, 'defaultAxesBox', 'on')

% Tiled layout 
set(groot,'defaultTiledLayoutPadding', 'tight')
set(groot,'defaultTiledLayoutTileSpacing', 'tight')

% Axis fontsize

FontSize = 20;
set(groot, 'defaultAxesFontSize', FontSize);
set(groot, 'defaultLegendFontSize', FontSize);
set(groot, 'defaultColorbarFontSize', FontSize);
set(groot, 'defaultTextarrowshapeFontSize', FontSize);
set(groot, 'defaultTextboxshapeFontSize', FontSize);
set(groot, 'defaultConstantLineFontSize', FontSize);
set(groot, 'defaultTextFontSize', FontSize);

%% other
% axis line thickness
set(groot, 'DefaultAxesLineWidth', 1);
set(groot, 'DefaultLegendAutoUpdate', 'off');

% legend object line thickness
set(groot, 'DefaultLineLineWidth', 3.0);
set(groot, 'DefaultLineMarkerSize', 15);
set(groot, 'DefaultStairLineWidth', 3.0);
set(groot, 'DefaultStairMarkerSize', 15);

% constant line (xline, yline) color
set(groot, 'defaultConstantLineColor', 'k');
set(groot, 'defaultConstantLineLineWidth', 3.0);
% default color map
% cmap = spring(128);
% set(groot, 'defaultFigureColormap', cmap);

% plot color default
% corder = [0,0,1;1,0,0;0.8,0.8,0;0.66015625,0.66015625,0.66015625;0,0,0;1,0.64453125,0;1,0,1;0,0.5,0.5;0,0,0.54296875;0,0.390625,0;0,1,1;0.59765625,0.1953125,0.796875];
% set(groot, 'defaultAxesColorOrder', corder);

% figure size
% set(groot, 'defaultFigureUnits','pixels')
% set(groot, 'defaultFigurePosition',[100 100 650 400])
% set(groot, 'defaultFigurePosition',[100 100 400 400])
width = 20; %cm
hwratio = 1;
set(groot, 'defaultFigureUnits', 'centimeters')
set(groot, 'defaultFigurePosition', [3 3 width, hwratio*width])

% grid on
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')

% legend settings
% set(groot, 'defaultLegendBox', 'on')
% set(groot, 'defaultLegendItemTokenSize', [10, 10])

% GUI (figure) settings
% set(groot, 'defaultTextboxshapeEdgeColor', [1 1 1])
set(groot, 'defaultTextboxshapeLineStyle', 'none')

% axes clipping (looks better than the default '3dbox')
% set(groot, 'defaultAxesClippingStyle', 'rectangle')

% set(groot, 'DefaultAxesXLimitmethod', 'tight')

set(groot,'DefaultFigureWindowStyle','docked')
% set(groot,'DefaultFigureWindowStyle','normal')

% 
% use this command to get all graphics objects
% types = unique(get(findall(gcf, '-property', 'Type'), 'Type'));

%% Notes on saving plots
%{
Want to make final fontsize and linewidth the same regardless of figure size
To make plots consistent with matlab appearance and saved png, use export_fig
JGCD: 2-column for final publication but can accomodate full row 

- subfigure with 2 columns:
    - LineWidth = 5, fontsize = 15


%}
