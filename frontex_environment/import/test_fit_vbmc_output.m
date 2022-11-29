% test_fit_output
clear all

% INPUT ------------------
modeltype = 'noisyKF_cfrule';
samplename = 'sample2';
load(sprintf('./sample_in/%s/out_fit_%s_ALL.mat',samplename,modeltype));
addpath(genpath('../../toolbox/fit_functions/vbmc-master'))
npar = 4; % some are 3 some are 4, specify explicitly here
ncond = 2;
% -----------------
nsubj = size(out_vbmc,1);


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
elseif strcmpi(modeltype,'noisyKF_cfrule')
    parstr = {'alpha','delta','zeta','tau'};
end


%% plot parameters
addpath ../../toolbox/plot_functions/
load ../../constants/condrgb.mat % loads condrgb

for ipar = 1:npar
    props.XLabel = parstr{ipar};
    plotViolins(squeeze(pars(idx,ipar,:)),condrgb*255,props,condstr);
end

%% plot parameter correlations
x = [pars(idx,:,1) pars(idx,:,2)];
corrtype = 'Spearman';
[r,p] = corr(x,'Type',corrtype);
dat.R = r;
dat.P = p;
dat.CorrType = corrtype;
if strcmpi(modeltype,'noisyKF')
    dat.labels = {'\alpha_B','\zeta_B','\tau_B','\alpha_F','\zeta_F','\tau_F'};
elseif strcmpi(modeltype,'noisyINF')
    dat.labels = {'h_B','\sigma_{I_B}','\sigma_{S_B}','h_F','\sigma_{I_F}','\sigma_{S_F}'};
elseif strcmpi(modeltype,'noisyKF_cfrule')
    dat.labels = {'\alpha_B','\delta_B','\zeta_B','\tau_B','\alpha_F','\delta_F','\zeta_F','\tau_F'};
end
plotCorrMatrix(dat);

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

