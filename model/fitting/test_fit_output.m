% test_fit_output
clear all

modeltype = 'noisyKF';
load(sprintf('../../frontex_environment/import/sample_in/pilot07/out_fit_%s_ALL.mat',modeltype));
addpath(genpath('../toolbox/vbmc-master'))

nsubj = size(out_vbmc,1);
npar = 3;
ncond = 2;

pars = nan(nsubj,npar,ncond);
for isubj = 1:nsubj
    if isempty(out_vbmc{isubj,1})
        continue
    end
    for icond = 1:2
        pars(isubj,:,icond) = out_vbmc{isubj,icond}.xavg;
    end
end

idx = ~isnan(pars(:,1,1));
condstr = {'bandit','fairy'};
if strcmpi(modeltype,'noisyKF')
    parstr = {'alpha','zeta','tau'};
elseif strcmpi(modeltype,'noisyINF')
    parstr = {'h','sigma','tau'};
end

%% plot parameters
addpath ../../toolbox/plot_functions/
load ../../constants/condrgb.mat % loads condrgb

ipar = 3;
props.XLabel = parstr{ipar};
plotViolins(squeeze(pars(idx,ipar,:)),condrgb*255,props,condstr)


%% plot parameter correlations
x = [pars(idx,:,1) pars(idx,:,2)];
[r,p] = corr(x);
figure
imagesc(r)
hold on
colormap(parula);
cb = colorbar;
clim([-1 1]);
if strcmpi(modeltype,'noisyKF')
    tklabels = {'\alpha_B','\zeta_B','\tau_B','\alpha_F','\zeta_F','\tau_F'};
elseif strcmpi(modeltype,'noisyINF')
    tklabels = {'h_B','\sigma_B','\tau_B','h_F','\sigma_F','\tau_F'};
end
xticklabels(tklabels);
yticklabels(tklabels);

%% examine posterior distribution of parameters for each subject/condition
isubj = 3;
icond = 1;
for isubj = 1:nsubj 
    if isempty(out_vbmc{isubj,1})
        continue
    end
    close all

    vbmc_plot(out_vbmc{isubj,icond}.vp);
    pause(.5)
end

