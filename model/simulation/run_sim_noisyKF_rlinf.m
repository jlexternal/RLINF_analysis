 % run_sim_noisyKF_rlinf
clear all

% set to true to check that the model produces identical output between 
% both conditions when noiseless
checkDeterministicCondition = false;

% -- Input: Specify source information --
samplename  = 'sample2';
modelkernel = 'noisyKF';
npar        = 3;
nsim        = 10;
savedir     = 'out'; % blank if current dir
savename    = sprintf('out_sim_noisyKF_nsim%d_%s',nsim,samplename);
%
%---------------------------------------------

% load fit parameters (location of parameter fits/vbmc output)
fitdir  = sprintf('../../frontex_environment/import/sample_in/%s',samplename); 
fitname = sprintf('out_fit_%s_ALL',modelkernel);
load(sprintf('%s/%s.mat',fitdir,fitname));

% load experiment data
load(sprintf('../../processed/%s/preprocessed_data_%s.mat',samplename,samplename));

nsubj = size(out_vbmc,1);
ncond = 2;

% extract parameters from VBMC output
pars = nan(nsubj,npar,ncond);
for isubj = 1:nsubj
    if isempty(out_vbmc{isubj,1})
        continue
    end
    for icond = 1:2
        pars(isubj,:,icond) = out_vbmc{isubj,icond}.xavg;
    end
end

if checkDeterministicCondition
    pars(:,1,2) = pars(:,1,1);
    pars(:,2:3,:) = 0;
end

idx = ~isnan(pars(:,1,1));
condstr = {'bandit','fairy'};
parstr = {'alpha','zeta','tau'};
nstype = 'weber';

addpath(genpath('../../toolbox/fit_functions'));

out_sim = cell(nsubj,ncond);
% run simulation
fprintf('Running simulations...\n');
for isubj = 1:nsubj
    if mod(isubj,10) == 0
        fprintf('Simulating KF model on subject %d...\n',isubj);
    end
    if isempty(out_vbmc{isubj,1})
        continue
    end
    cfg = [];
    cfg.nsim = nsim;

    % single out task vectors for specified subject
    dat = struct;
    dat.cond    = idx_cond(isubj,:);
    dat.blk     = idx_blk(isubj,:);
    dat.epi     = idx_epi(isubj,:);
    dat.trl     = idx_trial(isubj,:);
    dat.resp    = idx_blmn(isubj,:);
    
    % organize vectors for blue/moon feedback values
    bmst    = idx_bmstate(isubj,:);
    fbcorr  = idx_fbabs(isubj,:)/100;
    fb_ci   = [fbcorr' 1-fbcorr'];
    % value of MOON / BLUE color strength
    fb1     = fb_ci(sub2ind(size(fb_ci),(1:size(fb_ci,1))',-bmst'+2));
    
    for icond = 1:2
        cfg.condstr = condstr{icond};
        cfg.alpha   = pars(isubj,1,icond);
        cfg.zeta    = pars(isubj,2,icond);
        cfg.nstype  = nstype;
        cfg.tau     = pars(isubj,3,icond);

        % get relevant subject data
        trl = dat.trl(dat.cond == icond-1);
        blk = dat.blk(dat.cond == icond-1); 
        % order trials within block
        nt = sum(blk == 1);
        nb = max(blk);
        trl = trl - nt*(blk-1);

        idx_resp1 = trl == icond; % 1 for bandit, 2 for fairy
        resp1 = dat.resp(dat.cond == icond-1);
        resp1 = resp1(idx_resp1);
        resp1(resp1 == 0) = 2;
    
        cfg.r1      = fb1(dat.cond == icond-1);
        cfg.trl     = trl'; 
        cfg.firstresp = resp1;

        out_sim{isubj,icond} = sim_noisyKF_rlinf(cfg);

        bmstate = bmst(dat.cond == icond-1);
        bmstate(bmstate == 0) = 2;
        out_sim{isubj,icond}.corr = repmat(bmstate,[nsim 1]) == out_sim{isubj,icond}.rt;
%         fprintf('Accuracy cond %d: %.2f\n',icond,mean(out_sim{isubj,icond}.corr,'all'));
    end
end

% save structure
sim_out = struct;
sim_out.samplename  = samplename;
sim_out.nsim        = nsim;
sim_out.pars        = pars;
sim_out.fxname      = sprintf('sim_%s_rlinf',modelkernel);
sim_out.date        = datetime; 
sim_out.out         = out_sim;

% save to file
fullsavedir = sprintf('./%s/%s.mat',savedir,savename);
if isfile(fullsavedir)
    warning('File with name ''%s'' already exists!\n',savedir,savename);
    fprintf('Continue? (press any key to continue | Ctrl+C/Stop to terminate.\n');
    pause;
end
save(fullsavedir,'sim_out','-v7.3');
fprintf('File saved!\n');

%% deterministic model check

rts = nan(nsubj,292,2);
for isubj = 1:nsubj
    if isempty(sim_out.out{isubj,1})
        continue
    end
    for icond = 1:2
        pts(isubj,:,icond) = sim_out.out{isubj,icond}.pt;
    end
end

idx = ~isnan(pts(:,1,1));
test = pts(idx,:,1)-pts(idx,:,2);
teststr = {'passed','failed'};
fprintf('Deterministic test results: %s!\n',teststr{sum(test,'all')==0});