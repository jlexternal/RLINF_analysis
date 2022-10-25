function plotCorrMatrix(dat)
% Inputs:
%   dat: configuration file with following fields
%       R:          (matrix)    values of correlation coefficients
%  opt) P:          (matrix)    p-values for each correlation
%  opt) labels:     (cell arr)  labels of input values 
%  opt) CMap:       (string)    colormap (default: parula)
%  opt) CLim:       (array)     2x1 array for colorbar limits
%  opt) printflag:  (bool)      true if saving to disk 
%  opt) savename:   (string)    name of output file
%  opt) savedir:    (string)    relative directory of output file
%  opt) figprops:   (struct)    properties for final figure
%  opt) istriang:   (bool)      true if lower triangular matrix
%
% Jun Seok Lee - October 2022

if nargin < 1
    error('Missing data configuration!');
end
if ~isfield(dat,'R')
    error('Correlation matrix values are missing! (field ''R'')');
end
R = dat.R;
if any(isnan(R))
    error('Warning: NaNs found in input dataset. Proceeding...')
end
if isfield(dat,'P')
    P = dat.P';
    if any(P>1,'all') | any(P<0,'all')
        error('P-values found outside measure limits!');
    end
end
istriang = false;
if isfield(dat,'istriang')
    istriang = true;
end
npar = size(R,1);
if isfield(dat,'labels')
    labels = dat.labels;
end
cmap = 'parula';
if isfield(dat,'CMap')
    cmap = dat.CMap;
end
cLim = [-1 1];
if isfield(dat,'CLim')
    cLim = dat.CLim;
end
pbar = 1;
if isfield(dat,'figprops')
    figprops = dat.figprops;
    if isfield(figprops,'BoxAspectRatio'); pbar = figprops.BoxAspectRatio; end
    if isfield(figprops,'XLabel'); xlabelstr = figprops.XLabel; end
    if isfield(figprops,'YLabel'); ylabelstr = figprops.YLabel; end
end

printflag = false;
if isfield(dat,'printflag')
    printflag = true;
end

fprintf('Producing correlation matrix...\n');

if printflag
    if ~isfield(dat,'savename') | ~isfield(dat,'savedir')
        error('Name of file or target directory not provided (fields: ''savename'' or ''savedir'')');
    end
    hf = figure('Color','white');
else
    figure;
end
clf
% plot
imagesc(R);
hold on
for ipar = 1:npar
    for jpar = 1:npar
        if istriang & jpar <= ipar
            continue
        end
        text(ipar,jpar,sprintf('%s',sigstars(P(ipar,jpar))),...
            'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14);
    end
end
% plot styling
colormap(cmap);
cb = colorbar;
clim(cLim);

% figure styling
if exist('labels','var')
    xticks(1:npar);
    yticks(1:npar);
    xticklabels(labels);
    yticklabels(labels);
end
if exist('xlabelstr','var')
    xlabel(xlabelstr);
end
if exist('ylabelstr','var')
    ylabel(ylabelstr);
end
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);

if printflag
    set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
    set(gca,'FontName','Helvetica','FontSize',7.2);
    
    axes = findobj(gcf, 'type', 'axes');
    for a = 1:length(axes)
        if axes(a).YColor <= [1 1 1]
            axes(a).YColor = [0 0 0];
        end
        if axes(a).XColor <= [1 1 1]
            axes(a).XColor = [0 0 0];
        end
    end

    if ~isempty(figh)
        % save figure to file
        set(hf,'PaperPositionMode','manual', ...
            'PaperPosition',[2.5,13,16,figh],'PaperUnits','centimeters', ...
            'PaperType','A4','PaperOrientation','portrait');
        figure(hf);
        print(sprintf('%s/%s',savedir,savename),'-painters','-dpdf');
    end
end
end

function starstr = sigstars(pval)
    if pval < 0.001
        starstr = '***';
    elseif pval < 0.01
        starstr = '**';
    elseif pval < 0.05
        starstr = '*';
    else
        starstr = 'n.s.';
    end
end