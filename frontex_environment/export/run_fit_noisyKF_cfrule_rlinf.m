function run_fit_noisyKF_cfrule_rlinf(isubj)

% User defined input----------------------
samplename = 'sample1';
load(sprintf('./data/%s/preprocessed_data_%s.mat',samplename,samplename));
savekernel = 'out_fit_noisyKF_cfrule';
% ----------------------------------------

% load optimizers
addpath(genpath('./toolbox/fit_functions/'));

out_bads = cell(1,2);
out_vbmc = cell(1,2);

condstr = {'bandit','fairy'};

if ~isnan(idx_blmn(isubj,1))
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
    fb1     = fb_ci(sub2ind(size(fb_ci),(1:size(fb_ci,1))',-bmst'+2));
    dat.rt  = [fb1 1-fb1]; % value of MOON / BLUE color strength

    cfg = struct;
    for icond = 0:1
        cfg.condstr = condstr{icond+1};

        % get relevant subject data
        trl = dat.trl(dat.cond == icond);
        blk = dat.blk(dat.cond == icond); 

        % order trials within block
        nt = sum(blk == 1);
        nb = max(blk);
        trl = trl - nt*(blk-1);
        
        cfg.trl = trl';
        cfg.resp    = dat.resp(dat.cond == icond);
        cfg.rt      = dat.rt(dat.cond == icond,:);
        cfg.resp(cfg.resp == 0) = 2;
    
        cfg.fitalgo = 'bads';   % BADS algorithm
        cfg.nrun    = 1;       % 10 random starting points
        cfg.nsmp    = 1e3;      % 1e3 samples used by the particle filter
        cfg.nres    = 1e2;      % 1e2 validation samples used by BADS
    
        cfg.cfrule  = false;    % counterfactual rule flag (true or false)
        cfg.nstype  = 'weber';  % noise type (weber or white)
        cfg.chrule  = 'softm';  % choice rule (thomp or softm)
        %cfg.delta   = 0;       % comment out when not forcing delta to 0
    
        cfg.verbose = 2;        % 0, 1(basic), 2(super verbose) 
        cfg.noprior = true;     % for BADS only
    
        % fit BADS
        fprintf('Fitting subject %d over condition %d (BADS)...\n',isubj,icond);
        out_bads{isubj,icond+1} = fit_noisyKF_cfrule_rlinf(cfg);
    
    
        cfg.fitalgo = 'vbmc';
        cfg.pini.alpha = out_bads{isubj,icond+1}.alpha;
        cfg.pini.delta = out_bads{isubj,icond+1}.delta;
        cfg.pini.zeta  = out_bads{isubj,icond+1}.zeta;
        cfg.pini.tau   = out_bads{isubj,icond+1}.tau;
        cfg.noprior    = false;
    
        % fit VBMC
        fprintf('Fitting subject %d over condition %d (VBMC)...\n',isubj,icond);
        out_vbmc{isubj,icond+1} = fit_noisyKF_cfrule_rlinf(cfg);
    end
    
    dirroot = 'res';
    dirname = sprintf('./%s/%s',dirroot,samplename);

    % save file
    save(sprintf('%s/%s_subj%03d.mat',dirname,savekernel,isubj),'out_bads','out_vbmc');

end

end % function end