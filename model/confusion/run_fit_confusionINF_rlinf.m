function run_fit_confusionINF_rlinf(iset)

addpath(genpath('./toolbox/fit_functions/'));

samplename = 'pilot07';
load(sprintf('./data/%s/sim_dat_4confusion.mat',samplename),'sim');

dat = sim.dat(iset,:); % 2 element set array

out_bads = cell(1,2);
out_vbmc = cell(1,2);

condstr = {'bandit','fairy'};
for icond = 1:2
    cfg = struct;
    % get task data from sim structure
    cfg.condstr = condstr{icond};
    cfg.trl = dat{icond}.cfg.trl;
    cfg.blk = dat{icond}.cfg.blk;
    cfg.rt1 = dat{icond}.cfg.r1;
    % get simulated responses
    cfg.resp = dat{icond}.rt;

    cfg.fitalgo = 'bads';   % BADS algorithm
    cfg.nrun    = 10;        % 10 random starting points
    cfg.nsmp    = 1e3;      % 1e3 samples used by the particle filter
    cfg.nres    = 1e2;      % 1e2 validation samples used by BADS

    cfg.nstype  = dat{icond}.cfg.nstype;  % noise type (weber or white)
    cfg.verbose = 2;        % 0, 1(basic), 2(super verbose) 
    cfg.noprior = true;     % for BADS only

    % fit BADS
    fprintf('Fitting parameter set %d over the %s task (BADS)...\n',iset,condstr{icond});
    out_bads{iset,icond} = fit_noisyINF_rlinf(cfg);

    cfg.fitalgo = 'vbmc';
    cfg.pini.h          = out_bads{iset,icond}.h;
    cfg.pini.sigma_inf  = out_bads{iset,icond}.sigma_inf;
    cfg.pini.sigma_sel  = out_bads{iset,icond}.sigma_sel;
    cfg.noprior    = false;

    % fit VBMC
    fprintf('Fitting subject %d over the %s task (VBMC)...\n',iset,condstr{icond});
    out_vbmc{iset,icond} = fit_noisyINF_rlinf(cfg);
end

% save output
dirroot = 'res';
dirname = sprintf('./%s/%s',dirroot,samplename);
if not(isfolder(dirname))
    fprintf('Making folder %s in directory %s...\n',samplename,dirroot);
end

out_fit = struct;
out_fit.out_bads    = out_bads;
out_fit.out_vbmc    = out_vbmc;
out_fit.samplename  = samplename;
save(sprintf('%s/out_fit_noisyINF_conf_s%03d.mat',dirname,iset),'out_fit');

end