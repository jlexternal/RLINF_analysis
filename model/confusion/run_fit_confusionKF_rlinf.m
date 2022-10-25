function run_fit_confusionKF_rlinf(iset,lesionstr)

addpath(genpath('./toolbox/fit_functions/'));
islesion = false;
if exist('lesionstr','var')
    if ~ischar(lesionstr)
        error('Provide a string value for the lesioned parameter in 2nd argument!')
    elseif ismember(lesionstr,{'zeta','tau'})
        islesion = true;
    else
        error('Provided lesioned parameter not found in model!');
    end
end
samplename = 'pilot07';

% load appropriate data file (full sim, lesioned zeta OR tau)
if ~islesion
    load(sprintf('./data/%s/sim_dat_4confusion_KF.mat',samplename),'sim');
    fprintf('Loading data from full simulations...\n');
else
    load(sprintf('./data/%s/sim_dat_4confusion_no%s_KF.mat',samplename,lesionstr),'sim');
    fprintf('Loading data from lesioned %s simulations...\n',lesionstr);
end
lesionedpars = {'zeta','tau','none'};

dat = sim.dat(iset,:); % 2 element set array
    
for irun = 1:3
    out_bads = cell(1,3);
    out_vbmc = cell(1,3);
    
    condstr = {'bandit','fairy'};
    for icond = 1:2
        cfg = struct;
        % get task data from sim structure
        cfg.condstr = condstr{icond};
        cfg.trl = dat{icond}.cfg.trl;
        cfg.blk = dat{icond}.cfg.blk;
        cfg.rt1 = dat{icond}.cfg.r1/100;
        % get simulated responses
        cfg.resp = dat{icond}.rt;
    
        cfg.fitalgo = 'bads';   % BADS algorithm
        cfg.nrun    = 10;        % 10 random starting points
        cfg.nsmp    = 1e3;      % 1e3 samples used by the particle filter
        cfg.nres    = 1e2;      % 1e2 validation samples used by BADS
    
        cfg.nstype  = dat{icond}.cfg.nstype;  % noise type (weber or white)
        cfg.verbose = 2;        % 0, 1(basic), 2(super verbose) 
        cfg.noprior = true;     % for BADS only

        % fixed parameters for lesions
        if irun == 1 % lesion zeta
            cfg.zeta = 0;
            fprintf('Fitting model with no zeta...\n');
        elseif irun == 2
            cfg.tau = 0;
            fprintf('Fitting model with no tau...\n');
        else
            fprintf('Fitting full model...\n');
        end
    
        % fit BADS
        fprintf('Fitting task set %d over the %s task (BADS)...\n',iset,condstr{icond});
        out_bads = fit_noisyKF_rlinf(cfg);
    
        cfg.fitalgo = 'vbmc';
        cfg.pini.alpha = out_bads.alpha;
        cfg.pini.zeta  = out_bads.zeta;
        cfg.pini.tau   = out_bads.tau;
        cfg.noprior    = false;
    
        % fit VBMC
        fprintf('Fitting subject %d over the %s task (VBMC)...\n',iset,condstr{icond});
        out_vbmc = fit_noisyKF_rlinf(cfg);
    end
    
    % save output
    dirroot = 'res';
    dirname = sprintf('./%s/%s',dirroot,samplename);
    if not(isfolder(dirname))
        fprintf('Making folder %s in directory %s...\n',samplename,dirroot);
    end
    
    out_fit = struct;
    out_fit.lesionedpar = lesionedpars{irun};
    out_fit.out_bads    = out_bads;
    out_fit.out_vbmc    = out_vbmc;
    out_fit.samplename  = samplename;
    if ~islesion
        save(sprintf('%s/out_fit_noisyKF_conf_fullData_%sLesionedModel_s%03d.mat',dirname,lesionedpars{irun},iset),'out_fit');
    else
	    save(sprintf('%s/out_fit_noisyKF_conf_no%sData_%sLesionedModel_s%03d.mat',dirname,lesionstr,lesionedpars{irun},iset),'out_fit');
    end
end

end
