% compare_f
clear all

%  ----------- INPUT -----------
samplename  = 'sample2';
kernelarr   = {'out_fit_noisyKF_ALL',...
               'out_fit_noisyKF_cfrule_ALL'};
momenttype   = 'xavg'; % xavg or xmap
nsubj = 247; % to be replaced from input from constants LATER
% -------------------------------

condstr = {'bandit','fairy'};
nmod = numel(kernelarr);

load(sprintf('../../constants/constants_rlinf_%s.mat',samplename),'ncnd'); % load constants
load(sprintf('../../processed/%s/preprocessed_data_%s.mat',samplename,samplename),'idx_blmn');
load(sprintf('../../processed/%s/idx_TaskQuesAll.mat',samplename));

elbos = nan(nsubj,ncnd,nmod);
pars = struct;

% extract ELBO values from each model
for imod = 1:nmod
    load(sprintf('./sample_in/%s/%s.mat',samplename,kernelarr{imod}));
    
    pars.(kernelarr{imod}) = struct;
    for isubj = 1:nsubj
        if isempty(out_vbmc{isubj,1})
            continue
        end

        for icond = 1:ncnd
            if ~isfield(pars.(kernelarr{imod}),'pars')
                npars = numel(out_vbmc{isubj,icond}.(momenttype));
                pars.(kernelarr{imod}).pars = nan(nsubj,npars,ncnd);
            end
            elbos(isubj,icond,imod) = out_vbmc{isubj,icond}.elbo; % get ELBO
            pars.(kernelarr{imod}).pars(isubj,:,icond) = out_vbmc{isubj,icond}.(momenttype); % get parameters
        end
    end
end

%% Fixed effects BMS

model_rank = nan(nsubj,nmod,ncnd); 
for icond = 1:ncnd
    for isubj = 1:nsubj
        if isnan(elbos(isubj,1,1)) | idx_fullAll(isubj) == 0
            continue
        end
        for icond = 1:ncnd
            [~,model_rank(isubj,:,icond)] = sort(squeeze(elbos(isubj,icond,:)),'descend');
        end
    end
end

% model frequency
n_full = sum(idx_fullAll);
fprintf('\nModel frequencies: \n');
for imod = 1:nmod
    fprintf('%s\n',kernelarr{imod});
    for icond = 1:2
        fprintf('%s: %.2f  ',condstr{icond},(sum(model_rank(:,imod,icond) == 1))/n_full);
    end
    fprintf('\n\n');
end

for icond = 1:ncnd
    fprintf('Δ(ELBO)_%s: %.2f\n',condstr{icond},sum(elbos(idx_fullAll,icond,1)) - sum(elbos(idx_fullAll,icond,2)));
end

%% Parameter correlation matrix
addpath ../../toolbox/plot_functions/

for imod = 1:nmod
    par_mod = pars.(kernelarr{imod}).pars;
    npars = size(par_mod,2);

    for icond = 1:2
        x = [par_mod(idx_fullAll,:,1) par_mod(idx_fullAll,:,2)];
    end

    corrtype = 'Pearson';
    [r,p] = corr(x,'Type',corrtype);
    dat.R = r;
    dat.P = p;
    dat.CorrType = corrtype;
    % for KF model only
    if npars == 3
        dat.labels = {'\alpha_B','\zeta_B','\tau_B','\alpha_F','\zeta_F','\tau_F'};
    elseif npars == 4
        dat.labels = {'\alpha_B','\delta_B','\zeta_B','\tau_B','\alpha_F','\delta_F','\zeta_F','\tau_F'};
    end
    plotCorrMatrix(dat);
end





