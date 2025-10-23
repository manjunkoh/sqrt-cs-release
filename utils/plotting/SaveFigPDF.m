%% save pdf
% input >  fig : figure handle , filename : Filename before .pdf
%
% Example
% f1 = figure;
% plot(1:10);
% SaveFigPDF2(f1,'test')
%
% 
function y = SaveFigPDF(fig,filename,option,saveFileType)
if iscell(fig)
	fig = fig{1};
end
if nargin > 2
	if iscell(option)
		option = option{1};
	end
	if isfield(option, 'tight')
		if option.tight == 1
			set(fig.CurrentAxes, 'LooseInset', get(fig.CurrentAxes, 'TightInset'));
		else
			set(fig.CurrentAxes, 'LooseInset', [0.13, 0.11, 0.095, 0.075]);
		end
	else
		set(fig.CurrentAxes, 'LooseInset', get(fig.CurrentAxes, 'TightInset'));
	end
	if isfield(option, 'Flag_3D')
		Flag_3D = option.Flag_3D;
	else
		Flag_3D = 0;
	end
	if isfield(option, 'BackgroundColor')
		BackgroundColor = option.BackgroundColor;
	else
		BackgroundColor = [1, 1, 1]; % default
	end
	if isfield(option, 'ContentType')
		ContentType = option.ContentType;
	else
		ContentType = 'auto'; % default
	end
else
	set(fig.CurrentAxes, 'LooseInset', get(fig.CurrentAxes, 'TightInset'));
	Flag_3D = 0;
	BackgroundColor = [1, 1, 1]; % default
end
if exist('saveFileType', 'var')
	if ~iscell(saveFileType)
		saveFileType = {saveFileType};
    end
elseif isfield(option, 'saveFileType')
	if ~iscell(option.saveFileType)
		saveFileType = {option.saveFileType};
    else
        saveFileType = option.saveFileType;
    end
else
	saveFileType = {'all'};
end

set(fig,'Units','Inches');
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

for ii=1:numel(saveFileType)
	if strcmp(saveFileType{ii},'all') || strcmp(saveFileType{ii}, 'fig')
		saveas(fig, [filename, '.fig'])
	end
	if Flag_3D==0
		if strcmp(saveFileType{ii},'all') || strcmp(saveFileType{ii}, 'pdf')
			exportgraphics(fig,[filename, '.pdf'], 'ContentType', ContentType, 'BackgroundColor', BackgroundColor)
		end
		if strcmp(saveFileType{ii},'all') || strcmp(saveFileType{ii}, 'png')
			exportgraphics(fig,[filename, '.png'], 'BackgroundColor', BackgroundColor)
		end
	elseif Flag_3D==1
		if isfield(option, 'skewAngle')
			skewAngle = option.skewAngle;
		elseif isfield(option, 'angle')
			skewAngle = option.angle;
		else
			skewAngle = [-50, 8];
		end
		if numel(skewAngle)==2
			view(fig.CurrentAxes, skewAngle(1), skewAngle(2));
		elseif numel(skewAngle)==3
			view(fig.CurrentAxes, skewAngle);
		else
			error('invalid form: skewAngle')
		end
		set(fig.CurrentAxes, 'LooseInset', get(fig.CurrentAxes, 'TightInset'));
		if strcmp(saveFileType{ii},'all') || strcmp(saveFileType{ii}, 'png')
			exportgraphics(fig, [filename, '_angle', '.png'], 'Resolution', 300, 'BackgroundColor', BackgroundColor)
		end
		if strcmp(saveFileType{ii},'all') || strcmp(saveFileType{ii}, 'pdf')
			exportgraphics(fig, [filename, '_angle', '.pdf'], 'Resolution', 300, 'BackgroundColor', BackgroundColor)
		end
		view(fig.CurrentAxes, 0, 0);
		set(fig.CurrentAxes, 'LooseInset', get(fig.CurrentAxes, 'TightInset'));
		if strcmp(saveFileType{ii},'all') || strcmp(saveFileType{ii}, 'png')
			exportgraphics(fig, [filename, '_xz', '.png'], 'Resolution', 300, 'BackgroundColor', BackgroundColor)
		end
		if strcmp(saveFileType{ii},'all') || strcmp(saveFileType{ii}, 'pdf')
			exportgraphics(fig, [filename, '_xz', '.pdf'], 'ContentType', ContentType, 'BackgroundColor', BackgroundColor)
		end
		view(fig.CurrentAxes, 0, 90);
		set(fig.CurrentAxes, 'LooseInset', get(fig.CurrentAxes, 'TightInset'));
		if strcmp(saveFileType{ii},'all') || strcmp(saveFileType{ii}, 'png')
			exportgraphics(fig, [filename, '_xy', '.png'], 'Resolution', 300, 'BackgroundColor', BackgroundColor)
		end
		if strcmp(saveFileType{ii},'all') || strcmp(saveFileType{ii}, 'pdf')
			exportgraphics(fig, [filename, '_xy', '.pdf'], 'ContentType', ContentType, 'BackgroundColor', BackgroundColor)
		end
		view(fig.CurrentAxes, 90, 0);
		set(fig.CurrentAxes, 'LooseInset', get(fig.CurrentAxes, 'TightInset'));
		if strcmp(saveFileType{ii},'all') || strcmp(saveFileType{ii}, 'png')
			exportgraphics(fig, [filename, '_yz', '.png'], 'Resolution', 300, 'BackgroundColor', BackgroundColor)
		end
		if strcmp(saveFileType{ii},'all') || strcmp(saveFileType{ii}, 'pdf')
			exportgraphics(fig, [filename, '_yz', '.pdf'], 'ContentType', ContentType, 'BackgroundColor', BackgroundColor)
		end
	end
end % end if
end