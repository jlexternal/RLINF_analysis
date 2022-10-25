function run_fit_noisyKF_rlinf(isubj)

% User defined input----------------------
samplename = 'pilot07';
load(sprintf('./data/%s/preprocessed_data_%s.mat',samplename,samplename));
savekernel = 'out_fit_noisyKF';
% ------------------------------------

is_sim_recovery = true;
if is_sim_recovery
    load(sprintf('./data/%s/out_sim_4recov.mat',samplename));
end

addpath(genpath('./toolbox/fit_functions/'));

out_bads = cell(1,2);
out_vbmc = cell(1,2);

condstr = {'bandit', 'fairy'};

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
    dat.rt1 = fb1; % value of MOON / BLUE color strength
    
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
    
        
        if ~is_sim_recovery 
            % subject responses
            cfg.resp    = dat.resp(dat.cond == icond)';
        else
            sim_resps   = sim_out.out{isubj,icond+1}.rt';
            sim_resps(sim_resps==2) = 0;
            cfg.resp    = sim_resps;
        end
        cfg.trl     = trl';
        cfg.rt1     = dat.rt1(dat.cond == icond);
    
        cfg.fitalgo = 'bads';   % BADS algorithm
        cfg.nrun    = 1;        % 10 random starting points
        cfg.nsmp    = 1e3;      % 1e3 samples used by the particle filter
        cfg.nres    = 1e2;      % 1e2 validation samples used by BADS
    
        cfg.nstype  = 'weber';  % noise type (weber or white)
        cfg.verbose = 2;        % 0, 1(basic), 2(super verbose) 
        cfg.noprior = true;     % for BADS only
        
        % force learning or choice policy
        if false 
            % cfg.zeta = 0; % exact learning
            % cfg.tau = 0;  % argmax choice policy
        end
    
        % fit BADS
        fprintf('Fitting subject %d over the %s task (BADS)...\n',isubj,condstr{icond+1});
        out_bads{isubj,icond+1} = fit_noisyKF_rlinf(cfg);
    
        cfg.fitalgo = 'vbmc';
        cfg.pini.alpha = out_bads{isubj,icond+1}.alpha;
        cfg.pini.zeta  = out_bads{isubj,icond+1}.zeta;
        cfg.pini.tau   = out_bads{isubj,icond+1}.tau;
        cfg.noprior    = false;
    
        % fit VBMC
        fprintf('Fitting subject %d over the %s task (VBMC)...\n',isubj,condstr{icond+1});
        out_vbmc{isubj,icond+1} = fit_noisyKF_rlinf(cfg);
    end

    % save output
    dirroot = 'res';
    dirname = sprintf('./%s/%s',dirroot,samplename);
    if not(isfolder(dirname))
        fprintf('Making folder %s in directory %s...\n',samplename,dirroot);
    end
    
    if ~is_sim_recovery
        out_fit = struct;
        out_fit.out_bads    = out_bads;
        out_fit.out_vbmc    = out_vbmc;
        out_fit.samplename  = samplename;
        save(sprintf('%s/%s_s%03d.mat',dirname,savekernel,isubj),'out_fit');
    else
        out_rec = struct;
        out_rec.out_bads    = out_bads;
        out_rec.out_vbmc    = out_vbmc;
        out_rec.samplename  = samplename;
        save(sprintf('%s/%s_recov_s%03d.mat',dirname,savekernel,isubj),'out_rec');
    end
end

end % function end
