% organize_confusion
clear all

addpath ../simulation/ % add directory containing sim functions
nset = 100; % number of parameter sets to sample
imodel = 1; % 1-KF, 2-INF
nstype = 'weber'; % applicable for KF model only

%% 1) Generate 100 tasks
%       > From gen_trials_test.m (n=100; get 'fbs' and 'idx_epi' and save to file)
load('./tasks/tasks_n100_4confusion.mat','idx_epi','fbs');

% organize task data to fit simulation code
% Requires: cfg
%              .trl/ trial index within block
%              .blk/ block index within task
%              .r1/  value of option 1 (blue/moon)
epiPerBlk = 6;
ntrPerBlk = 73; % extra trial added
idx_blk = [];
for iblk = 1:4
    idx_blk = [idx_blk iblk*ones(1,ntrPerBlk)];
end
idx_fb1 = fbs; clearvars fbs;
idx_blk = repmat(idx_blk,[nset 1]);
idx_trl = repmat(repmat(1:ntrPerBlk,[1 4]),[nset 1]);    

%% 2) Sample 100 sets of parameters from the subjects' fitted parameters
              %   1          2
modeltypes = {'noisyKF','noisyINF'};
modelpars = {'alpha','zeta','tau'; 'h','sigma_inf','sigma_sel'};
load(sprintf('../../processed/pilot07/out_fit_%s.mat',modeltypes{imodel}),'out');

pars_subj = out.pars;
idx_subj = ~isnan(pars_subj(:,1,1));
npar   = size(pars_subj,2);
parset = nan(nset,npar);

% 2a) Organize sets of non-correlated parameters
isCorrelated = true;
while isCorrelated
    for ipar = 1:npar
        % fit distributions
        if ipar == 1
            distname = 'Beta';
            pd = fitdist(reshape(pars_subj(idx_subj,ipar,:),[numel(pars_subj(idx_subj,ipar,:)) 1]),distname);
        else
            distname = 'Normal';
            pd = fitdist(reshape(log(pars_subj(idx_subj,ipar,:)),[numel(pars_subj(idx_subj,ipar,:)) 1]),distname);
        end
        if ipar == 1 
            parset(:,ipar) = betarnd(pd.a,pd.b,nset,1);
        else
            parset(:,ipar) = exp(normrnd(pd.mu,pd.sigma,nset,1));
        end
    end

    % check for correlation and repeat if found
    % > Should give null correlation matrix (correlate set to itself, diagonal meaningless)
    [~,p] = corr(parset);
    isCorrelated = any(p < 0.05,'all');
end

% visualize parameter distribution of simulations and subjects
clf
for ipar = 1:3
    subplot(3,1,ipar);
    histogram(parset(:,ipar),'Normalization','pdf');
    hold on
    histogram(reshape(pars_subj(idx_subj,ipar,:),[numel(pars_subj(idx_subj,ipar,:)) 1]),'Normalization','pdf');
end

%% 3) Simulate behavior from the 100 above. 
nsim = 1;

% lesioned simulations
islesioned = true;
pars_fixed = [.4 1 0;...    % no lesions
              .4 0 0.05;... % lesioned zeta
              .4 1 0];      % lesioned tau

irun = 3; % choose lesion set

%             1        2
condstr = {'bandit','fairy'};

out_sim = cell(nset,2);
fprintf('Running simulations using %s...\n',modeltypes{imodel});
for iset = 1:nset
        fprintf('Running simulations on set number %d...\n',iset);
    for icond = 1:2
        r1 = [];
        for ib = 1:4
            if icond == 1 % bandit : last trial NaN
                r1  = [r1 idx_fb1(iset,1+(72*(ib-1)):ib*72) NaN];
            else % fairy  : first trial NaN
                r1 = [r1 NaN idx_fb1(iset,1+(72*(ib-1)):ib*72)];
            end
        end
        if imodel == 1
            r1 = r1/100;
        end

        cfg = struct;
        for ipar = 1:3
            cfg.(sprintf('%s',modelpars{imodel,ipar})) = pars_fixed(irun,ipar);
        end
        
        cfg.condstr = condstr{icond};
        cfg.nstype  = nstype;
        cfg.nsim    = nsim;
        cfg.r1      = r1;
        cfg.blk     = idx_blk(iset,:);
        cfg.trl     = idx_trl(iset,:);
        
        if strcmpi(modeltypes{imodel},'noisyKF')
            out_sim{iset,icond} = sim_noisyKF_rlinf(cfg);
        elseif strcmpi(modeltypes{imodel},'noisyINF')
            out_sim{iset,icond} = sim_noisyINF_rlinf(cfg);
        else
            error('Unexpected value of modeltype!');
        end
        out_sim{iset,icond}.cfg = cfg;
    end
end

sim = struct;
sim.dat = out_sim;
sim.genpars = parset;
sim.description = 'Simulated data for confusion matrix with constant alpha, zeta; lesioned tau';
sim.modeltype = modeltypes{imodel};
fprintf('Description for simulation is as follows:\n %s\n\n',sim.description);
fprintf('Continue with save?\n\n');
pause

savename = 'sim_dat_4confusion_notau_KF';
save(savename,'sim');
fprintf('Simulation data saved in %s.mat !\n',savename);

%% 4) Fit the model on simulated data (in server)

%%%%%%%%%%%%

%% 5) Produce confusion matrix 
%       > Should give correlation matrix with only the diagonal 
clear all
% load fits and organize
% correlate with generative parameter sets

% Inputs:
imodel = 1; % 1-KF, 2-INF
nset = 100; % match from above

modeltypes = {'KF','INF'};
momenttype = 'xavg';
modelpars = {'alpha','zeta','tau'; 'h','sigma_inf','sigma_sel'};

load(sprintf('sim_dat_4confusion_%s.mat',modeltypes{imodel}));

pars_sim = nan(nset,3,2);
pars_rec = nan(nset,3,2);
for isubj = 1:nset
    
    load(sprintf('./fits/out_fit_noisy%s_conf_s%03d.mat',modeltypes{imodel},isubj));

    out_vbmc(isubj,:) = out_fit.out_vbmc(isubj,:);

    for icond = 1:2
        pars_rec(isubj,:,icond) = out_vbmc{isubj,icond}.(momenttype);
        for ipar = 1:3
            pars_sim(isubj,ipar,icond) = sim.dat{isubj,icond}.(modelpars{imodel,ipar});
        end
    end
end

