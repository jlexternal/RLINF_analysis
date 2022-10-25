function run_fit_noisyKF_rlinf(isubj)

samplename = 'pilot07';
load(sprintf('../../processed/%s/preprocessed_data_%s.mat',samplename,samplename));

is_sim_recovery = false;
if is_sim_recovery
    load ./out_sim_4recov.mat
end

addpath(genpath('../toolbox/fit_functions/'));

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
            sim_resps   = sim_out.out{isubj,icond}.rt';
            sim_resps(sim_resps==2) = 0;
            cfg.resp    = sim_resps;
        end
        cfg.trl     = trl';
        cfg.rt1     = dat.rt1(dat.cond == icond);

        cfg.fitalgo = 'bads';   % BADS algorithm
        cfg.nrun    = 1;        % 10 random starting points
        cfg.nsmp    = 1e3;      % 1e3 samples used by the particle filter
        cfg.nres    = 1e2;      % 1e2 validation samples used by BADS
    
        cfg.verbose = 2;        % 0, 1(basic), 2(super verbose) 
        cfg.noprior = true;     % for BADS only
        
        % force learning or choice policy
        if false 
            % cfg.sigma_inf = 0; % exact inference
            % cfg.sigma_sel = 0; % argmax choice policy
        end
    
        % fit BADS
        fprintf('Fitting subject %d over the %s task (BADS)...\n',isubj,condstr{icond+1});
        out_bads{isubj,icond+1} = fit_noisyKF_rlinf(cfg);
    
        cfg.fitalgo = 'vbmc';
        cfg.pini.alpha     = out_bads{isubj,icond+1}.alpha;
        cfg.pini.sigma_inf = out_bads{isubj,icond+1}.sigma_inf;
        cfg.pini.sigma_sel = out_bads{isubj,icond+1}.sigma_sel;
        cfg.noprior    = false;
    
        % fit VBMC
        fprintf('Fitting subject %d over the %s task (VBMC)...\n',isubj,condstr{icond+1});
        out_vbmc{isubj,icond+1} = fit_noisyKF_rlinf(cfg);
    end
    
    if ~is_sim_recovery
        out_fit = struct;
        out_fit.out_bads = out_bads;
        out_fit.out_vbmc = out_vbmc;
        out_fit.samplename = samplename;
        save('out_fit_test.mat','out_fit');
    else
        out_rec = struct;
        out_rec.out_bads = out_bads;
        out_rec.out_vbmc = out_vbmc;
        out_rec.samplename = samplename;
        save('out_fit_recov.mat','out_rec');
    end
end

end % function end