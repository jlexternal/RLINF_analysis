% run_mode_tester
%
% Simulates some stereotypical data generated on the simulation function
% corresponding to some model. 
%
% 1/ Failure to meet expected outcome of simulation represents misspecification 
% or improper model.
% 2/ Failure to fit represents either irrecoverability of the model or an error 
% in either the simulation or fit function.
%
% Jun Seok Lee - Nov 2022
%
clear all
clc
addpath ../toolbox/cprintf/ % colored console output text

% ------ Input: Model ID ----------------
modeltype   = 'noisyKF_cfrule';
projectname = 'rlinf';
% ---------------------------------------

% Add directories for simulation and fits and task
addpath(genpath('./simulation'));
addpath(genpath('./fitting'));
addpath(genpath('../toolbox/fit_functions'));
load('./reliability/task_testcases_rlinf.mat','testcases'); % loads stereotypical task data
ncase = numel(testcases);
fprintf('Task test data loaded (%d test cases).\n',ncase);

% get information about the specific model
modelkernel = sprintf('%s_%s',modeltype,projectname);
simfnstr    = sprintf('sim_%s',modelkernel);
fitfnstr    = sprintf('fit_%s',modelkernel);
simfn = str2func(simfnstr); % localize model simulation function
fitfn = str2func(fitfnstr); % localize model fitting function
fprintf('Loading information about model (%s)...\n\n',modeltype);
cfg = struct;
cfg.getInfo = true;
simfn(cfg); % get model information

% ---------------- Input: ----------------
% set stereotypical parameter sets of input model from info above
parstr = {'alpha','delta','zeta','tau'};
isparnoisy = logical([0 0 1 1]); % flag noise parameters
nsmp = 1e3; % number of samples for particle filter in fits
nsim = 100; % number of simulations if model is noisy
condstr = {'bandit','fairy'};
nstype = 'weber';
cfrule = false;
chrule = 'softm';
% ---------------------------------------
ncond = numel(condstr);

%% 1. Specify model and task data

parset = cell(0);

% ---------------- Input: ----------------
% parset{1} = struct;
% parset{1}.description  = 'No learning (alpha=0)';
% parset{1}.pars         = [0 0 0 0];

% parset{1} = struct;
% parset{1}.description  = 'Full learning (alpha=1)';
% parset{1}.pars         = [1 0 0 0];
% 
% parset{2} = struct;
% parset{2}.description  = 'Medium learning (alpha=0.5); no decay (delta=0)';
% parset{2}.pars         = [.5 0 0 0];
% 
% parset{3} = struct;
% parset{3}.description  = 'Medium learning (alpha=0.5); full decay (delta=1)';
% parset{3}.pars         = [.5 1 0 0];
% 
parset{1} = struct;
parset{1}.description  = 'Medium learning, standard decay, standard noise (alpha=0.5, delta=0.2, zeta=1.0)';
parset{1}.pars         = [.5 .2 1 0];

parset{2} = struct;
parset{2}.description  = 'Medium learning, standard decay, softmax policy (alpha=0.5, delta=0.2, tau=0.1)';
parset{2}.pars         = [.5 .2 0 .1];

parset{3} = struct;
parset{3}.description  = 'Medium learning, standard decay, noise, softmax policy (alpha=0.5, delta=0.2, zeta=1.0 tau=0.1)';
parset{3}.pars         = [.5 .2 1 .1];

nparset = numel(parset);
% ----------------------------------------

% run simulation and check sim output for test success
out_sim = cell(nparset,ncase);

for iparset = 1:nparset
    fprintf('Simulating parameters (%s) on test cases...\n',parset{iparset}.description);
    % set up simulation config
    cfg = struct;
    checkDeterministicCondition = true; % check between conditions

    % load model parameter set
    cfg.nsim = 1;
    for ipar = 1:numel(parstr)
        cfg.(parstr{ipar}) = parset{iparset}.pars(ipar);
        if isparnoisy(ipar) && parset{iparset}.pars(ipar) ~= 0
            cfg.nsim = nsim; % multiple simulations for non-deterministic models
            checkDeterministicCondition = false;
        end
    end
    % load model settings
    cfg.nstype = nstype;
    cfg.cfrule = cfrule;

    % load task data
    for icase = 1:ncase
        cfg.r1  = testcases{icase}.idx_fb; % reward value of option 1
        cfg.trl = testcases{icase}.idx_trl; % trial index
        nb = sum(cfg.trl == 1);
        cfg.firstresp = ones(nb,1); % fix first response to 1

        for icond = 1:ncond
            cfg.condstr = condstr{icond};
            out_sim{iparset,icase}{icond} = simfn(cfg);
            out_sim{iparset,icase}{icond}.r1 = cfg.r1;
            out_sim{iparset,icase}{icond}.trl = cfg.trl;
            out_sim{iparset,icase}{icond}.checkDeterministicCondition = checkDeterministicCondition;
        end
    end
end
cprintf('green','All simulations done.\n');

%% 2. Check simulation output

% Input: choose simulations to analyze
set2analyze = 3% 1:3; %1:nparset;
% -------------------------------------

fprintf('\n');
fprintf('Press ''P'' or ''F'' to indicate visual test Pass/Fail..\n');
for iparset = set2analyze
    fprintf('Checking output of simulated responses from parameter specification:\n');
    fprintf('--------------------------------------------------------------------\n');
    cprintf('*text',' %s\n',parset{iparset}.description);
    fprintf('--------------------------------------------------------------------\n');
    
    % go through the test cases
    for icase = 1:ncase
        cprintf('*text','Test case: ');
        cprintf('_text','%s\n',testcases{icase}.description);
        
        trl = out_sim{iparset,icase}{1}.trl;
        nt_all = numel(trl);
        % ---------------- Input: ----------------
        % Note: These variables depend on the model. Specify the variables
        %       required by the analysis.
        nsim = size(out_sim{iparset,icase}{icond}.rt,1);
        rt   = nan(nsim,nt_all,2);   % nsim, trial, condition
        pt   = nan(nsim,nt_all,2);
        r_ch = nan(nsim,nt_all,2);
        mt   = nan(nsim,nt_all,2,2); % nsim, trial, option, condition
        for icond = 1:ncond
            rt(:,:,icond)   = out_sim{iparset,icase}{icond}.rt; % responses from condition
            pt(:,:,icond)   = out_sim{iparset,icase}{icond}.pt; % response 1 probability
            r_ch(:,:,icond) = out_sim{iparset,icase}{icond}.r_ch; % feedback received
            mt(:,:,:,icond) = out_sim{iparset,icase}{icond}.mt; % tracked values
        end
        % ----------------------------------------
        r_ch(rt == 2)   = 1-r_ch(rt == 2); % sign the feedbacks for visualization
        rt(rt == 2) = 0; % convert 2 to 0 for visualization purposes below
        
        % check if the two conditions produce identical results if model deterministic
        detTestStr = '  Deterministic model test between conditions: ';
        if out_sim{iparset,icase}{1}.checkDeterministicCondition
            idx_rt_1 = true(1,nt_all);
            idx_rt_2 = true(1,nt_all);

            idx_excl_1 = [find(trl == 1)-1 nt_all]; % excluded trials for bandit
            idx_excl_1(idx_excl_1==0) = [];
            idx_excl_2 = find(trl == 1); % excluded trials for fairy
            idx_rt_1(idx_excl_1) = false;
            idx_rt_2(idx_excl_2) = false;

            % index matching across conditions: rt_1(1:end-1) == rt_2(2:end)
            fprintf('%s',detTestStr);
            if all(rt(1,idx_rt_1,1) == rt(1,idx_rt_2,2))
                cprintf('green','Passed!\n');
            else
                cprintf('red','Failed!\n');
            end
            clearvars idx_rt_1 idx_rt_2 idx_excl_1 idx_excl_2
        end

        % check trajectory of simulations with the task values
        visTestStr = '  Visual model check: ';
        fprintf(pad(visTestStr,numel(detTestStr),'right'));

        figure(1);
        clf
        sgtitle(sprintf('Model simulation check\nModel spec: %s\nTask spec: %s',...
                parset{iparset}.description,testcases{icase}.description));
        for icond = 1:ncond
            % feedback seen, tracked values, option chosen
            subplot(ncond,1,icond);
            title(sprintf('Condition: %s',condstr{icond}),'FontSize',12);
            for it = find(trl == 1)
                xline(it-.5,'HandleVisibility','off');
            end
            ylim([-.1 1.1]);
            xlim([1 nt_all+4]);
            yline(.5,':','HandleVisibility','off');
            hold on
            yline(0,'Color','green','LineWidth',2,'Alpha',.5,'HandleVisibility','off');
            yline(1,'Color','red','LineWidth',2,'Alpha',.5,'HandleVisibility','off');

            % visualize means
            scatter([1:nt_all]-.5,mean(r_ch(:,:,icond),1),'x');           % feedback
            scatter([1:nt_all]-.1,mean(mt(:,:,1,icond),1),'o','red');     % tracked option 1
            scatter([1:nt_all]+.1,mean(mt(:,:,2,icond),1),'o','green');   % tracked option 2
            scatter([1:nt_all],mean(pt(:,:,icond),1),100,'square');       % probability of option 1
            scatter(1:nt_all,mean(rt(:,:,icond),1),'ko','filled');        % response given
            % legend
            legend({'Feedback','Tracker 1','Tracker 2','p(option1)','Choice'},'Location','southeast');
        end
        while true
            isKeyPress = waitforbuttonpress;
            keyValue = double(get(gcf,'CurrentCharacter'));
            if ~isKeyPress
                continue
            elseif ismember(keyValue,[112 80]) % p or P
                cprintf('green','Passed!\n');
                break
            elseif ismember(keyValue,[102 70]) % f or F
                cprintf('red','Failed!\n');
                break
            elseif ismember(keyValue,[117 85])
                cprintf('yellow','Unclear\n');
                break
            end 
        end
    end
    fprintf('\n');
end

%% 3. Run model fitter on simulated data for fit success

% Input: Fitting code parameters
nrun = 10;  % initial starting points for BADS algorithm
verbose = 0; % 0-'none', 1-'final', 2-'iter' (for bads/vbmc) 
parsets = 1:nparset;
cases = 1:ncase;
% -------------------------------------

if exist('out_vbmc') == 1
    warning('Variable name ''out_vbmc'' already exists in the workspace!');
    warning('Override? (press any key to continue)');
    pause
end

% data structures for fit output
out_bads = cell(nparset,ncase,ncond);
out_vbmc = cell(nparset,ncase,ncond);

% loop through the different sets of parameter values previously simulated
parfor iparset = parsets
    % load task data and simulated responses
    for icase = cases
        % set input structure for fit
        cfg = struct;
        r1  = testcases{icase}.idx_fb';
        cfg.rt = [r1 1-r1];
        cfg.trl = testcases{icase}.idx_trl;

        % load model settings
        cfg.nstype = nstype;
        cfg.cfrule = cfrule;
        cfg.chrule = chrule;
        cfg.verbose = verbose;

        % fit conditions
        for icond = 1:ncond
            cfg.condstr = condstr{icond};
            % load simulated responses
            out_sim_struct = out_sim{iparset,icase}{icond};
            % load only one of the simulated response trajectories if noisy
            itraj = 1;
            if ~out_sim_struct.checkDeterministicCondition
                itraj = randi(nsim);
            end
            cfg.resp = out_sim_struct.rt(itraj,:);
            
            % run BADS fitting algorithm
            cfg.fitalgo = 'bads';
            cfg.nsmp    = nsmp;
            cfg.nrun    = nrun;
            cfg.noprior = true;

            fprintf('Running BADS algorithm on %s over %s...\n',parset{iparset}.description,testcases{icase}.description);
            out_bads{iparset,icase,icond} = fitfn(cfg);

            % run VBMC fitting algorithm
            cfg.fitalgo = 'vbmc';
            % load point-estimate values as prior in VBMC
            for ipar = 1:numel(parstr)
                cfg.pini.(parstr{ipar}) = out_bads{iparset,icase,icond}.(parstr{ipar});
            end
            cfg.noprior = false;
            out_vbmc{iparset,icase,icond} = fitfn(cfg);

            fprintf('Fitting parameters \n  (%s) over test case\n   ',parset{iparset}.description);
            fprintf('%s\n COMPLETE\n',testcases{icase}.description);
        end
    end
end

% save the output
save(sprintf('./reliability/out_fits_%s.mat',string(datetime('now','Format','ddMMyyyy_HHmm'))),'out_bads','out_vbmc');

%% 4. Check parameter recovery
% Input: Choose generative parameter set to examine
iparset = 3;
% ------------------------

ylims = [0 1; 0 1; 0 3; 0 .2];
% check parameters
pars = parset{iparset}.pars;
fprintf('\n');
figure(2);
clf
sgtitle(sprintf(['%s\n' ...
    'line: target, x:xmap, o:xavg+xstd'],parset{iparset}.description));
% get fit parameters
for icase = 1:ncase
    fprintf('%s\n',testcases{icase}.description);
    xdisp = [-.1 .1];
    for icond = 1:2
        xmap_fit = out_vbmc{iparset,icase,icond}.xmap;
        xavg_fit = out_vbmc{iparset,icase,icond}.xavg;
        xstd_fit = out_vbmc{iparset,icase,icond}.xstd;
        
        for ipar = 1:numel(parstr)
            subplot(ncase,numel(parstr),(icase-1)*4+ipar);
            hold on
            scatter(1+xdisp(icond),xavg_fit(ipar),'bo');
            scatter(1+xdisp(icond),xmap_fit(ipar),'rx');
            plot((1+xdisp(icond))*ones(2,1),xavg_fit(ipar)+([1 -1]*xstd_fit(ipar)),'Color','b');

            % cosmetics
            if icond == 2
                elbos = [out_vbmc{iparset,icase,1}.elbo out_vbmc{iparset,icase,2}.elbo];
                plot(xdisp+1,pars(ipar)*ones(2,1),'k');
                xlim(1+xdisp*3);
                xticks(1+xdisp);
                xticklabels(condstr);
                ylim(ylims(ipar,:)+xdisp);
                if ipar == 1
                    ylabel(sprintf('ELBO:\n%.2f\n%.2f',elbos));
                    set(get(gca,'YLabel'),'Rotation',0,'HorizontalAlignment','right','VerticalAlignment','middle');
                end
                if icase == 1
                    title(sprintf('%s',parstr{ipar}));
                end
            end
        end
    end
end

%% 5. Check particle trajectories
% Input: ----------
iparset = 3;
icase = 4;
icond = 2;
momentstr = 'xavg';
sz_max = 40;
% -----------------
disp('');
fprintf('\nRunning particle filter on fits of %s\n',parset{iparset}.description);

% prepare particle filter config and output
out_pf = cell(ncase,ncond);
cfg = struct;
r1  = testcases{icase}.idx_fb';
cfg.rt = [r1 1-r1];
cfg.trl = testcases{icase}.idx_trl;
% load model settings
cfg.nstype = nstype;
cfg.cfrule = cfrule;
cfg.chrule = chrule;
cfg.nsmp   = nsmp;
% split by condition
cfg.condstr = condstr{icond};
% load simulated responses inputted into fit
cfg.resp = out_vbmc{iparset,icase,icond}.cfg.resp;
% set parameters
pars = out_vbmc{iparset,icase,icond}.(momentstr);
xnam = out_vbmc{iparset,icase,icond}.xnam;
for ipar = 1:numel(xnam)
    cfg.(xnam{ipar}) = pars(ipar);
end
cfg.noprior = false;
% run particle filter
out_pf{icase,icond} = fitfn(cfg);

% plot trajectories
figure(3);
clf
fitparstr = strjoin(strcat(parstr,cellstr(string(round(pars,2)))));
title(sprintf('%s\n%s\n%s\ncondition: %s',parset{iparset}.description,testcases{icase}.description,fitparstr,condstr{icond}));
mt = out_pf{icase,icond}.mt_raw;
wt = out_pf{icase,icond}.wt_raw;
rt = out_pf{icase,icond}.cfg.rt;
resp = cfg.resp;
nt = size(wt,1);
hold on
for it = 1:nt
    % size particles by weight
    wt_it = wt(it,:)/max(wt(it,:));
    sz = wt_it*sz_max;
    sz(sz==0) = 1e-3;
    % plot both options with x-offset
    scatter(it*ones(nsmp,1)-.1,squeeze(mt(it,1,:)),sz,'o','MarkerFaceColor','red','MarkerFaceAlpha',.05,'MarkerEdgeColor','none');
    scatter(it*ones(nsmp,1)+.1,squeeze(mt(it,2,:)),sz,'o','MarkerFaceColor','green','MarkerFaceAlpha',.05,'MarkerEdgeColor','none');
end
% trajectories of particles
plot((1:nt)-.1,squeeze(mt(:,1,:)),'Color',[1 0 0 0.01]);
plot((1:nt)+.1,squeeze(mt(:,2,:)),'Color',[0 1 0 0.01]);
% min/max of particles
scatter((1:nt)-.1,max(squeeze(mt(:,1,:)),[],2),100,'^','MarkerFaceColor','none','MarkerEdgeColor','red');
scatter((1:nt)-.1,min(squeeze(mt(:,1,:)),[],2),100,'v','MarkerFaceColor','none','MarkerEdgeColor','red');
scatter((1:nt)+.1,max(squeeze(mt(:,2,:)),[],2),100,'^','MarkerFaceColor','none','MarkerEdgeColor','green');
scatter((1:nt)+.1,min(squeeze(mt(:,2,:)),[],2),100,'v','MarkerFaceColor','none','MarkerEdgeColor','green');
for it = 1:nt
    % plot outcome/stimuli
    if icond == 1
        scatter((it-.5)+1,rt(it,resp(it)),80,'kx');
    elseif icond == 2
        scatter((it-.5)+1,rt(it,1),80,'kx');
    end
end
% plot responses
resp(resp==2)=0;
scatter(1:nt,resp,'ko','filled'); % response given

% cosmetics
yoffset = [-.1 .1];
ylimcnd = [0 1];
ymidcnd = 0.5;

ylim(ylimcnd+yoffset);
yline(ymidcnd,':');
yline(ylimcnd(1),'Color','green','LineWidth',2,'Alpha',.5,'HandleVisibility','off');
yline(ylimcnd(2),'Color','red','LineWidth',2,'Alpha',.5,'HandleVisibility','off');

