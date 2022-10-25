 % run_sim_noisyKF_rlinf
clear all

% set to true to check that the model produces identical output between 
% both conditions when noiseless
checkDeterministicCondition = false;

% specify fit dataset to simulate
samplename = 'pilot07';

% load fit parameters (location of parameter fits/vbmc output)
fitdir = sprintf('../../processed/%s',samplename);
fitname = 'out_fit_noisyINF';
load(sprintf('%s/%s.mat',fitdir,fitname));

nsim = 1;
load(sprintf('../../processed/%s/preprocessed_data_%s.mat',samplename,samplename));

% -- Input: Specify output file information --
savedir  = '../recovery/out'; % blank if current dir
% savename = sprintf('out_sim_noisyINF_nsim%d_%s',nsim,samplename);
savename = 'out_sim_noisyINF_4recov';
% --------------------------------------------

nsubj = size(out.pars,1);
npar = 3;
ncond = 2;

pars = out.pars;

if checkDeterministicCondition
    pars(:,1,2) = pars(:,1,1);
    pars(:,2:3,:) = 0;
end

idx = ~isnan(pars(:,1,1));
condstr = {'bandit','fairy'};
parstr = {'h','sigma_inf','sigma_sel'};

addpath(genpath('../../toolbox/fit_functions'));

out_sim = cell(nsubj,ncond);
% run simulation
fprintf('Running simulations...\n');
for isubj = 1:nsubj
    if mod(isubj,10) == 0
        fprintf('Simulating noisy Inference model on subject %d...\n',isubj);
    end
    if isnan(pars(isubj,1,1))
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
    fbcorr  = idx_fbabs(isubj,:);
    fb_ci   = [fbcorr' 100-fbcorr'];
    % value of MOON / BLUE color strength
    fb1     = fb_ci(sub2ind(size(fb_ci),(1:size(fb_ci,1))',-bmst'+2));
    
    for icond = 1:2
        cfg.condstr   = condstr{icond};
        cfg.h         = pars(isubj,1,icond);
        cfg.sigma_inf = pars(isubj,2,icond);
        cfg.sigma_sel = pars(isubj,3,icond);

        % get relevant subject data
        trl = dat.trl(dat.cond == icond-1);
        blk = dat.blk(dat.cond == icond-1); 
        % order trials within block
            nt = sum(blk == 1);
            nb = max(blk);
            trl = trl - nt*(blk-1);
    
        cfg.r1      = fb1(dat.cond == icond-1);
        cfg.trl     = trl';
        
        out_sim{isubj,icond} = sim_noisyINF_rlinf(cfg);
    end
end

% save structure
sim_out = struct;
sim_out.samplename  = samplename;
sim_out.nsim        = nsim;
sim_out.pars        = pars;
sim_out.fxname      = 'sim_noisyINF_rlinf';
sim_out.date        = datetime; 
sim_out.out         = out_sim;

% save to file
fullsavedir = sprintf('./%s/%s.mat',savedir,savename);
if isfile(fullsavedir)
    warning('File with name ''%s'' already exists!\n',savedir,savename);
    fprintf('Continue? (press any key to continue | Ctrl+C/Stop to terminate.\n');
    pause
end
save(fullsavedir,'sim_out');
fprintf('File saved!\n');

%% deterministic model check

rts = nan(23,292,2);
for isubj = 1:23
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
