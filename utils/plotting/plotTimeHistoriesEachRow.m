function [fig, figprop, ax] = plotTimeHistoriesEachRow(t_his, x_his, lineSpecs, labelList, fig, figprop, options)
%/*************************************************************************
%* File Name     : plotTimeHistoriesEachRow.m
%* Description   : description
%* Author        : Kenshiro (Ken) Oguri
%* Affiliation   : School of Aeronautics & Astronautics, Purdue University
%* Creation Date : Thu Sep  2 12:30:42 2021
%* Language      : Matlab
%* Version       : 1.0.0: first version
%                : 1.0.1: 
%**************************************************************************
%* INPUTS  :
%*      - input:
%*
%* OUTPUTS :
%*      - output:
%*
%*************************************************************************/

[Nx, Nt] = size(x_his);
if Nt ~= numel(t_his)
	error('t and u data not consistent')
end
if ~isempty(labelList)
	Nx = numel(labelList.y);
end

if ~exist('fig', 'var') || isempty(fig)
	isNewFig = 1;
	fig 	= figure();
	figprop = struct('tight', 0, 'Flag_3D', 0);
	ax = tiledlayout(Nx, 1, 'TileSpacing', 'tight');
else
	isNewFig = 0;
	figure(fig);
end
if isempty(lineSpecs)
	lineSpecs = struct('Type', 'plot', 'Line', '-', 'Color', 'k');
end
if ~exist('options', 'var')
    options = [];
end

% plot each component in a different row
for ii=1:Nx
	if isNewFig
		nexttile(ax)
	else
		nexttile(ii)
	end
	hold on
	if iscell(lineSpecs.Line)
		if isfield(lineSpecs, 'Type')
			switch lineSpecs.Type{ii}
		 	case 'stairs'
				p_u_i = stairs(t_his, x_his(ii,:), lineSpecs.Line{ii});
		 	case 'plot'
				p_u_i = plot(t_his, x_his(ii,:), lineSpecs.Line{ii});
			end
		else
			p_u_i = stairs(t_his, x_his(ii,:), lineSpecs.Line{ii});
		end
	else
		if isfield(lineSpecs, 'Type')
			switch lineSpecs.Type
		 	case 'stairs'
				p_u_i = stairs(t_his, x_his(ii,:), lineSpecs.Line);
		 	case 'plot'
				p_u_i = plot(t_his, x_his(ii,:), lineSpecs.Line);
			end
		else
			p_u_i = stairs(t_his, x_his(ii,:), lineSpecs.Line);
		end
    end
    if isfield(lineSpecs, 'Color')
    	if iscell(lineSpecs.Color)
	    	p_u_i.Color = lineSpecs.Color{ii};
    	else
	    	p_u_i.Color = lineSpecs.Color;
        end
    end
	if isfield(lineSpecs, 'LineWidth')
		if iscell(lineSpecs.LineWidth)
			p_u_i.LineWidth = lineSpecs.LineWidth{ii};
		else
			p_u_i.LineWidth = lineSpecs.LineWidth;
		end
	end
	if isfield(options, 'ylim')
		if iscell(options.ylim)
			if ~isempty(options.ylim{ii})
				ylim(options.ylim{ii});
			end
		else
			if ~isempty(options.ylim)
				ylim(options.ylim);
			end
		end
	end
	if isNewFig
		ylabel(labelList.y{ii}); grid on
		if ii==Nx
			xlabel(labelList.x)
		end
	end
end
% axis tight

% make axes same scale
xlim_ = [inf,-inf];
ylim_ = [inf,-inf];
try
for ii=1:Nx
	nexttile(ii)
	xlim_j = xlim;
	xlim_ = [min([xlim_(1), xlim_j(1)]), max([xlim_(2), xlim_j(2)])];
	ylim_j = ylim;
	ylim_ = [min([ylim_(1), ylim_j(1)]), max([ylim_(2), ylim_j(2)])];
end
end
if isfield(options, 'xlog') && options.xlog
	for ii=1:Nx
		nexttile(ii)
		set(gca, 'XScale', 'log')
	end
else
	% do nothing
end

if isfield(options, 'ylog') && options.ylog
	for ii=1:Nx
		nexttile(ii)
		set(gca, 'YScale', 'log')
	end
else
	% do nothing
end
try
if isfield(options, 'samexlim') && ~options.samexlim
	% do nothing
else
	for ii=1:Nx
		nexttile(ii)
		xlim(xlim_)
	end
end
end
try
if isfield(options, 'sameylim') && options.sameylim
	for ii=1:Nx
		nexttile(ii)
		ylim(ylim_)
	end
end
end
end