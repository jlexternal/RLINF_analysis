function plotViolins(X,rgbs,props,labels,printflag,savename,savedir)
% Inputs:
%   X:          data to plot, dim 1: source, dim 2: group
%   printflag:  true if exporting to file
%   rgbs:       array of RGB values (0-255) for each group
%   props:      structure with all the extra stuff
%   labels:     cell array of group labels
%   printflag:  true if saving to disk
%   savename:   name of output file
%   savedir:    relative directory of output file
%
% Jun Seok Lee - October 2022

if nargin < 2; error('Not enough arguments!'); end
    
% default values
mrksz       = 50;
datwid      = 0.2;
pbar        = .75;
momenttype  = 'mean';
errortype   = 'sem';
quantiles   = [.025 .975];

if exist('props')
    if isfield(props,'MarkerSize'); mrksz = props.MarkerSize; end
    if isfield(props,'ViolinWidth'); datwid = props.ViolinWidth; end
    if isfield(props,'BoxAspectRatio'); pbar = props.BoxAspectRatio; end
    if isfield(props,'MomentType'); momenttype = lower(props.MomentType); end % Mean/Median
    if isfield(props,'ErrorType'); errortype = lower(props.ErrorType); end % StDev/SEM/CI
    if strcmpi(errortype,'CI') & isfield(props,'Quantiles')
        quantiles   = props.Quantiles;  % 2 element vector of CI quantiles
    end
    if isfield(props,'XLabel'); xlabelstr = props.XLabel; end
    if isfield(props,'YLabel'); ylabelstr = props.YLabel; end
    if isfield(props,'YLim'); ylims = props.YLim; end
end

npts = size(X,1);
ngps = size(X,2);
if nargin < 4; labels = 1:ngps; end
if nargin < 5; printflag = false; end

linewidth = 2;
if printflag 
    if ~isfield(props,'FigureHeight')
        error('Height of figure not provided! (props.FigureHeight)');
    end
    figh = props.FigureHeight;
    if ~exist('savename')
        savename = sprintf('defaultViolinOut_%s',string(round(rand*10000)));
    end
    if ~exist('savedir')
        savedir = '.';
    end
    linewidth = 1;
end

if size(rgbs,1) ~= ngps
    error('Number of groups does not equal number of RGB sets! (%d != %d)',size(rgbs,1),ngps);
end

rgbs = rgbs/255;
if any(isnan(X),'all')
    fprintf('NaNs found in data! Continuing operation. Check input data array...\n');
end

fprintf('Producing violin plots with central moment (%s) and error (%s)...\n',momenttype, errortype);

if printflag
    hf = figure('Color','white');
else
    figure;
end
clf
hold on
for ig = 1:ngps
    x = X(:,ig);

    % violin
    xk = linspace(min(x),max(x),100);
    pk = ksdensity(x,xk);
    pk = pk/max(pk);
    str = interp1(xk,pk,x);
    jit = linspace(-1,+1,numel(x));
    for j = randperm(numel(x))
        p = scatter(ig+jit(j)*str(j)*datwid,x(j),mrksz*.8, ...
            'MarkerFaceColor',0.5*(rgbs(ig,:)+1),'MarkerEdgeColor','none');
        p.MarkerFaceAlpha = 0.4;
    end

    % error bars
    switch errortype
        case 'stdev'
            xerr = mean(x) + [-1 1]*std(x);
        case 'sem'
            xerr = mean(x) + [-1 1]*std(x)/sqrt(npts);
        case 'quantiles'
            xerr = quantile(x,quantiles);
        case 'quartiles'
            xerr = quantile(x,[.25 .75]);
    end
    plot([1 1]*ig,xerr,'Color',rgbs(ig,:),'LineWidth',linewidth);

    % moment 
    switch momenttype
        case 'mean'
            xmoment = mean(x);
        case 'median'
            xmoment = median(x);
    end
    scatter(ig,xmoment,mrksz,rgbs(ig,:),'LineWidth',linewidth,'MarkerFaceColor',[1 1 1]);
end
xticks(1:ngps);
xlim([.5 ngps+.5]);
xticklabels(labels);
if exist('ylims')
    ylim(ylims);
end
if exist('ylabelstr')
    ylabel(ylabelstr);
end
if exist('xlabelstr')
    xlabel(xlabelstr);
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