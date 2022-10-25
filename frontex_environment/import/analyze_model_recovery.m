% analyze_model_recovery
clear all

modeltype = 'INF';

addpath(genpath('../../toolbox/plot_functions'));

samplename = 'pilot07';
load(sprintf('../../constants/constants_rlinf_%s',samplename),'ncnd','condrgb'); % load  constants
condstr = {'bandit','fairy'};

% load original fit parameter
filename = 'pars_fit_noisyINF.mat';
load(sprintf('./sample_out/%s/%s',samplename,filename));
xpars = out.pars; % original parameters (x)
% parstr = {'alpha','zeta','tau'};
parstr = {'h','sigma_{inf}','sigma_{sel}'};
% prange = [0 .1 1; 0 .1 4; 0 .01 .4]; % min, increment, max of parameter ranges
prange = [0 .1 1; 0 .1 2; 0 .1 6]; % min, increment, max of parameter ranges
idx = ~isnan(xpars(:,1,1));

% load recovery
filename = 'pars_fit_noisyINF_recov.mat';
load(sprintf('./sample_out/%s/%s',samplename,filename));
ypars = out.pars; % comparison parameters (y)

fprintf('Comparing parameter fits (recovery)...\n');

npar = size(xpars,2);

%% check fit to recovered parameters

figure;
clf
ctr = 0;
for ipar = 1:npar
    fprintf('\nParameter recovery %s model (%s)...\n',modeltype,parstr{ipar});
    for icond = 1:ncnd
        ctr = ctr + 1;
        x = xpars(idx,ipar,icond);
        y = ypars(idx,ipar,icond);
        % stats
        [r,p] = corr(x,y,'type','Pearson');
        fprintf('r=%.4f, p=%.4f (%s)\n',r,p,condstr{icond});
        
        % plot
        xrange = prange(ipar,1):prange(ipar,2):prange(ipar,3);
        [pn,s] = polyfit(x,y,1);
        [py,d] = polyconf(pn,xrange,s,'alpha',0.05,'predopt','curve');
        subplot(npar,ncnd,ctr);
        hold on
        % confidence bands
        s = shadedErrorBar(xrange,py,d,'patchSaturation',.1,'lineprops',{'LineWidth',1,'Color',condrgb(icond,:)});
        set(s.edge,'LineStyle','none');
        % data
        scatter(x,y,20,'MarkerFaceColor',0.5*(condrgb(icond,:)+1),'MarkerEdgeColor','none');
        % identity
        plot(xrange,xrange,':');
        xlim(prange(ipar,[1 3]));
        ylim(prange(ipar,[1 3]));
        xlabel(sprintf('Original %s',parstr{ipar}));
        ylabel(sprintf('Recovered %s',parstr{ipar}));
        if ipar == 1
            title(sprintf('%s',condstr{icond}));
        end
    end
end
sgtitle(sprintf(['Parameter recoverability (model: %s-cf)\n' ...
    'Source: %s\nnsubj=%d'],modeltype,samplename,sum(idx)));

%% check inter-parameter correlations

parcombs = [1 2; 1 3; 2 3];
parlims = [0 0.8; 0 2; 0 6];
figure;
clf
ctr = 1;
for icomb = 1:3
    ipar = parcombs(icomb,1);
    jpar = parcombs(icomb,2);
    for icond = 1:ncnd
        subplot(3,2,ctr);
        hold on
        % original fits
        x = xpars(idx,ipar,icond);
        y = xpars(idx,jpar,icond);
        xrange = [min(x) max(x)];
        [pn,s] = polyfit(x,y,1);
        [py,d] = polyconf(pn,xrange,s,'alpha',0.05,'predopt','curve');
        s = shadedErrorBar(xrange,py,d,'patchSaturation',.1,'lineprops',{'LineWidth',1,'Color',condrgb(icond,:)});
        set(s.edge,'LineStyle','none');
        scatter(x,y,30,'MarkerFaceColor',0.5*(condrgb(icond,:)+1),'MarkerEdgeColor','none');
        xlim(parlims(ipar,:));
        ylim(parlims(jpar,:));
        [r1,p1] = corr(x,y,'type','Pearson');
    
        % recovered fit
        x = ypars(idx,ipar,icond);
        y = ypars(idx,jpar,icond);
        xrange = [min(x) max(x)];
        [pn,s] = polyfit(x,y,1);
        [py,d] = polyconf(pn,xrange,s,'alpha',0.05,'predopt','curve');
        s = shadedErrorBar(xrange,py,d,'patchSaturation',.1,'lineprops',{'--','LineWidth',1,'Color',[0 0 0]});
        set(s.edge,'LineStyle','none');
        scatter(x,y,30,'o','MarkerFaceColor',[0 0 0],'MarkerEdgeColor','none','MarkerFaceAlpha','.5');
        xlim(parlims(ipar,:));
        ylim(parlims(jpar,:));
        [r2,p2] = corr(x,y,'type','Pearson');

        xlabel(parstr{ipar});
        ylabel(parstr{jpar});
        title(sprintf('fit: r=%.4f, p=%.4f\nrec: r=%.4f, p=%.4f',r1,p1,r2,p2),'FontSize',12);
    
        ctr = ctr + 1;
    end
end
sgtitle(sprintf(['Interparameter correlation (model: %s-cf)\n' ...
    'Source: %s\nnsubj=%d'],modeltype,samplename,sum(idx)));

