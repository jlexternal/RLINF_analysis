% analyze_model_recovery
clear all

addpath(genpath('../toolbox/plot_functions'));

samplename = 'pilot07';
load(sprintf('../constants/constants_rlinf_%s',samplename),'ncnd','condrgb'); % load  constants

% load original fit parameter
filename = 'pars_fit_noisyINF';
load(sprintf('./fitting/out/%s/%s.mat',samplename,filename));

xpars = out.pars; % original parameters (x)
clearvars out
parstr = {'alpha','sigma_inf','sigma_sel'};
prange = [0 .1 1; 0 .1 4; 0 .01 .4]; % min, increment, max of parameter ranges
idx = ~isnan(xpars(:,1,1));

% load recovery
filename = 'pars_fit_rec_noisyINF';
load(sprintf('./fitting/out/%s/%s.mat',samplename,filename));
ypars = out.pars; % comparison parameters (y)
condstr = {'bandit','fairy'};
%%
fprintf('Comparing parameter fits (recovery)...\n');

npar = size(xpars,2);
% check fit to recovered parameter correlations
figure;
clf
ctr = 0;
for ipar = 1:npar
    fprintf('\nParameter recovery KF model (%s)...\n',parstr{ipar});
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
        xlabel(sprintf('Original %s',parstr{ipar}));
        ylabel(sprintf('Recovered %s',parstr{ipar}));
        if ipar == 1
            title(sprintf('%s',condstr{icond}));
        end
    end
end
sgtitle(sprintf(['Parameter recoverability (model: KF-cf)\n' ...
    'Source: %s\nnsubj=%d'],samplename,sum(idx)));

