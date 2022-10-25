function plotReversalCurves(dat)
% Inputs:
%   dat: configuration file with following fields
%       X1:         (matrix)    main data (iteration, trial, group)
%  opt) X2:         (optional)  comparison data
%  rec) revtrl:     (int)       trial index after reversal
%  rec) rgbs:       (matrix)    RGB values (0-255) for each group
%  rec) dlabel:     (cell)      array for dataset label (e.g. 'data','sim')
%  rec) glabel:     (cell)      array of group labels (e.g. 'cond1','cond2')
%  opt) printflag:  (bool)      true if saving to disk 
%  opt) savename:   (string)    name of output file
%  opt) savedir:    (string)    relative directory of output file
%  opt) patchsat:   (double)    opacity of error shading 
%  opt) lineprops:  (cell)      names and values of curve drawing properties 
%                               (e.g. {'LineWidth',2})
%  opt) figprops:   (struct)    properties for final figure
%
% Requirements:
%   shadedErrorBar by Rob Campbell (2009)
%
% Jun Seok Lee - October 2022

if nargin < 1
    error('Missing data configuration!');
end
if ~isfield(dat,'X1')
    error('Trial data is missing! (field ''X1'')');
end
X1 = dat.X1;
if any(isnan(X1))
    warning('Warning: NaNs found in input dataset. Proceeding...')
end
% comparasion dataset if desired
ndat = 1;
if isfield(dat,'X2')
    X2 = dat.X2;
    ndat = 2;
    if any(not(isnan(X1) == isnan(X2)))
        warning('Warning: Location(s) of NaNs in main data do not match of comparison!\n');
    end
end
% dimensions of data
[nsubj,ntrl,ngrp] = size(X1);
% location of trial after reversal
if isfield(dat,'revtrl')
    revtrl = dat.revtrl;
else
    revtrl = 1;
    fprintf('Reversal trial not provided. Plotting normal learning curves... (field ''revtrl'')');
end 
% get RGB color values
if isfield(dat,'rgbs')
    rgbs = dat.rgbs/255;
end

% defaults
patchsat    = 0.1;
linewidth   = 2;
pbar        = 1.25;

if ~isfield(dat,'lineprops')
    lineprops = {'LineWidth',linewidth};
else
    lineprops = dat.lineprops;
end
if isfield(dat,'patchsat')
    patchsat = dat.patchsat;
end

if ~isfield(dat,'dlabel')
    dlabel_default = {'Data1','Data2'};
    dlabel = dlabel_default{1:ndat};
else
    dlabel = dat.dlabel;
end
if ~isfield(dat,'glabel')
    glabel_default = {};
    for igrp = 1:ngrp
        glabel_default{end+1} = sprintf('Cond%02d',igrp);
    end
    glabel = glabel_default;
else
    glabel = dat.glabel;
end
if isfield(dat,'figprops')
    figprops = dat.figprops;
    if isfield(figprops,'BoxAspectRatio'); pbar = figprops.BoxAspectRatio; end
    if isfield(figprops,'XLabel'); xlabelstr = figprops.XLabel; end
    if isfield(figprops,'YLabel'); ylabelstr = figprops.YLabel; end
    if isfield(figprops,'YTicks'); yTicks = figprops.YTicks; end
    if isfield(figprops,'YLim'); ylims = figprops.YLim; end
end

printflag = false;
if isfield(dat,'printflag')
    printflag = true;
end

fprintf('Producing reversal/learning curves with curve (mean) and error (SEM)...\n');

if printflag
    if ~isfield(dat,'savename') | ~isfield(dat,'savedir')
        error('Name of file or target directory not provided (fields: ''savename'' or ''savedir'')');
    end
    hf = figure('Color','white');
else
    figure;
end
clf
hold on
% plot
for igrp = 1:ngrp
    cv    = nanmean(X1(:,:,igrp),1); % moment 
    cv_e  = nanstd(X1(:,:,igrp),[],1)/sqrt(nsubj); % error (SEM)

    % style curve
    lp_cv = lineprops;
    lp_cv{end+1} = 'Color'; lp_cv{end+1} = rgbs(igrp,:); % color
    
    % plot main data curve
    c1 = shadedErrorBar(1:ntrl,cv,cv_e,'patchSaturation',patchsat,'lineprops',lp_cv);
    set(c1.edge,'LineStyle','none');

    % plot comparison data (if provided)
    if ndat == 2
        cv_a    = nanmean(X2(:,:,igrp),1);
        cv_e_a  = nanstd(X2(:,:,igrp),[],1)/sqrt(nsubj);
        c2 = shadedErrorBar(1:ntrl,cv_a,cv_e_a,'patchSaturation',patchsat,'lineprops',['--',lp_cv]);
        set(c2.edge,'LineStyle','none');
    end
end
% figure styling
xticks(1:ntrl);
xticklabels([(1:revtrl-1)-revtrl 1:(ntrl-(revtrl-1))])
xline(revtrl-.5,':')
xlim([1 ntrl]);
legendstr = {};
for igrp = 1:ngrp
    for idat = 1:ndat
        legendstr{end+1} = sprintf('%s (%s)',glabel{igrp},dlabel{idat});
    end
end
legend(legendstr,'Location','southeast','FontSize',14);
if exist('ylims')
    ylim(ylims);
end
if exist('yticks')
    yticks(yTicks);
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