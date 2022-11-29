% behavioral_analysis
clear all
close all
clc

addpath ./toolbox/plot_functions/

% -------------------------INPUT--------------------------------
samplename = 'sample2'; % Pilot - 'pilot'
useFullAll = true;      % use subjects with full task and ques
comparesim = false;     % compare to simulation data 
comparert = true;

% for simulation analysis
modelkernel = 'noisyKF_cfrule';
nsim = 10;
simfilename = sprintf('processed_out_sim_%s_nsim%d_%s.mat',modelkernel,nsim,samplename);

% for retest analysis
if comparert
    dat_rt = load(sprintf('./processed/%s/preprocessed_data_retest_%s.mat',samplename,samplename)); % load the raw data structure for data sample
end
% --------------------------------------------------------------
if comparesim & comparert
    error('Choose either simulations or retest data to analyze (not both)!');
end

load(sprintf('./constants/constants_rlinf_%s.mat',samplename)); % load nblk, ncond, ntrl, samplename
load(sprintf('./processed/%s/preprocessed_data_%s.mat',samplename,samplename)); % load the raw data structure for data sample

nsubj = size(idx_fb,1);
nexcl = sum(isnan(idx_fb(:,1)));

idx_rt = ~isnan(dat_rt.idx_blmn(:,1));
idx_rt = cat(1,idx_rt,zeros(nsubj-numel(idx_rt),1));
idx_rt = logical(idx_rt);

fprintf('Number of subjects: %d\n',nsubj-nexcl);
condstr = {'bandit','fairy'};

idx_firsttrl = 1+73*(0:7);

% Compare participant behavior with simulations
if comparesim
    fprintf('\nSimulation comparison on (comparesim)\n');
    fprintf('Loading %s...\n',simfilename);
    load(sprintf('./processed/%s/%s',samplename,simfilename));
    sim_resps = out.resps;
    sim_resps(sim_resps == 2) = 0;
    nsim = out.nsim;
end

% Specify subjects included in analyses
if useFullAll
    load(sprintf('./processed/%s/idx_TaskQuesAll.mat',samplename),'idx_fullAll');
    idx_subj = idx_fullAll;
    if comparert
        load(sprintf('./processed/%s/idx_task_rt.mat',samplename));
        load(sprintf('./processed/%s/idx_ques_rt.mat',samplename));
        if size(idx_ques_rt) ~= size(idx_task_rt)
            idx_ques_rt = idx_ques_rt';
        end
        idx_subj = idx_task_rt & idx_ques_rt;
    end
else
    idx_subj = ~isnan(idx_fb(:,1));
end

% compute evidence factor
perr = 0.3;
x = -1:0.001:+1;
xfac = fzero(@(b)trapz(x(x > 0),1./(1+exp(x(x > 0)*b))/trapz(x,1./(1+exp(x*b))))-perr,[0,10]);
fprintf('evidence factor = %.3f\n',xfac);

% process feedback for simulation responses
if comparesim
    idx_fb_sim = nan(nsim,ntrl*ncnd*nblk,nsubj);
    for isubj = 1:nsubj
        if ~idx_subj(isubj)
            continue
        end
        resps = sim_resps(:,:,isubj); % blue/moon
        fbabs_mat = repmat(idx_fbabs(isubj,:),[nsim 1]);
        bmst_mat = repmat(idx_bmstate(isubj,:),[nsim 1]);

        corr = resps == bmst_mat; % convert to sim responses to isCorrect values
        fbabs_mat(~corr) = 100-fbabs_mat(~corr);
        idx_fb_sim(:,:,isubj) = fbabs_mat;
    end
    clearvars resps fbabs_mat bmst_mat corr
end

%% Plot predicted and observed outcome distribution
figure
hold on
xlim([-1 1]);
histogram(reshape(idx_fbabs(idx_fullAll,:),[],1)/100*2-1, ...
    'Normalization','pdf','NumBins',99,'FaceColor',[0.8,0.8,0.8],'EdgeColor','w');
plot(x,1./(1+exp(-x*xfac))/trapz(x,1./(1+exp(-x*xfac))),'r-','LineWidth',2);
hold off
xlabel('fbabs rescaled in [-1,+1]');
hold off

%% Overall accuracy and repetitions
pcorr = nan(nsubj,2);
prepe = nan(nsubj,2);

if comparert
    pcorr_rt = nan(nsubj,2);
    prepe_rt = nan(nsubj,2);
end

for isubj = 1:nsubj
    if ~idx_subj(isubj)
        continue
    end
    if comparert & ~idx_rt(isubj)
        continue
    end
    
    for icond = 0:1
        % Calculating accuracy: ignore the 1st trial of each block
        idx = idx_cond(isubj,:) == icond & ~ismember(idx_abstr(isubj,:),idx_firsttrl);
        pcorr(isubj,icond+1) = nanmean(idx_corr(isubj,idx));
        if comparert & idx_rt(isubj)
            idx = dat_rt.idx_cond(isubj,:) == icond & ~ismember(dat_rt.idx_abstr(isubj,:),idx_firsttrl);
            pcorr_rt(isubj,icond+1) = nanmean(dat_rt.idx_corr(isubj,idx));
        end
        
        % repetitions must be calculated within each block
        nrepe = 0;
        if comparert; nrepe_rt = 0; end
        for iblk = 1:4
            idx = idx_cond(isubj,:) == icond & idx_blk(isubj,:) == iblk;
            resps = idx_blmn(isubj,idx);
            if comparert
                idx = dat_rt.idx_cond(isubj,:) == icond & dat_rt.idx_blk(isubj,:) == iblk;
                resps_rt = dat_rt.idx_blmn(isubj,idx);
            end
            if icond == 0 % bandit (compare starting 2nd to 1st trial)
                nrepe = nrepe + sum(resps(2:end) == resps(1:end-1));
                if comparert; nrepe_rt = nrepe_rt + sum(resps_rt(2:end) == resps_rt(1:end-1)); end
            else % apples: (compare starting 3rd to 2nd trial, since first trial is nothing)
                nrepe = nrepe + sum(resps(3:end) == resps(2:end-1));
                if comparert; nrepe_rt = nrepe_rt + sum(resps_rt(3:end) == resps_rt(2:end-1)); end
            end
        end
        
        if icond == 0
            ntot = 72;
        else
            ntot = 71;
        end 
        prepe(isubj,icond+1) = nrepe/(ntot*4);
        if comparert; prepe_rt(isubj,icond+1) = nrepe_rt/(ntot*4); end
    end
end

if false
    save(sprintf('./processed/%s/pcor_%s.mat',samplename,samplename),'pcorr');
end

for icond = 1:2
    if comparert
        fprintf('General performance in retest(%s, N=%d)\n',samplename,sum(idx_rt));
        fprintf('Comparison of subjects who have done the retest...\n')
        fprintf('%s\n',condstr{icond});
        fprintf('Accuracy: %.2f ± %.2f\n',nanmean(pcorr(idx_rt,icond)),nanstd(pcorr(idx_rt,icond)));
        fprintf('Accuracy(retest): %.2f ± %.2f\n',nanmean(pcorr_rt(idx_rt,icond)),nanstd(pcorr_rt(idx_rt,icond)));

        fprintf('Repeats: %.2f ± %.2f\n',nanmean(prepe(idx_rt,icond)),nanstd(prepe(idx_rt,icond)));
        fprintf('Repeats(retest): %.2f ± %.2f\n\n',nanmean(prepe_rt(idx_rt,icond)),nanstd(prepe_rt(idx_rt,icond)));
        continue
    end
    fprintf('General performance (%s, N=%d)\n',samplename,nsubj-nexcl);
    fprintf('%s\n',condstr{icond});
    fprintf('Accuracy: %.2f ± %.2f\n',nanmean(pcorr(:,icond)),nanstd(pcorr(:,icond)));
    fprintf('Repeats: %.2f ± %.2f\n\n',nanmean(prepe(:,icond)),nanstd(prepe(:,icond)));
    
end

if ~comparert
    figure
    hold on
    title('Correlation p(switch)','FontSize',14)
    x = 1-prepe(:,1);
    y = 1-prepe(:,2);
    idx = ~isnan(x);
    [pn,s] = polyfit(x(idx),y(idx),1);
    xrange = 0:.01:1;
    [py,d] = polyconf(pn,xrange,s,'alpha',0.05,'predopt','curve');
    s = shadedErrorBar(xrange,py,d,'patchSaturation',.1,'lineprops',{'LineWidth',2});
    set(s.edge,'LineStyle','none');
    plot(0:1,0:1,':');
    scatter(x(idx),y(idx));
    xlim([0 .6]);
    ylim([0 .6]);
    xlabel('Bandit');
    ylabel('Fairy');
    [r,p] = corr(x(idx),y(idx));
    fprintf('Correlation p(switch): r = %.4f, p = %.4f\n',r,p);
    hold off
end

if comparert
    measstr = {'corr','switch'};
    figure
    clf
    sgtitle(sprintf('Retest\nN=%d',sum(idx_rt)));
    for imeas = 1:2
        for icond = 1:2
            fprintf('%s\n',condstr{icond});
            subplot(2,2,icond+(2*(imeas-1)));
            title(sprintf('p(%s)\n%s',measstr{imeas},condstr{icond}));
            hold on
            if imeas == 1
                x = pcorr(idx_rt,icond);
                y = pcorr_rt(idx_rt,icond);
            else
                x = 1-prepe(idx_rt,icond);
                y = 1-prepe_rt(idx_rt,icond);
            end
            scatter(x,y);
            [r,p] = corr(x,y);
            fprintf('Correlation p(%s): r = %.4f, p = %.4f\n',measstr{imeas},r,p);
            plot(0:1,0:1,':');
            [pn,s] = polyfit(x,y,1);
            xrange = min(x):.01:max(x);
            [py,d] = polyconf(pn,xrange,s,'alpha',0.05,'predopt','curve');
            s = shadedErrorBar(xrange,py,d,'patchSaturation',.1,'lineprops',{'LineWidth',2});
            set(s.edge,'LineStyle','none');
            xlim([.5 .9]-(imeas-1)*[.5 .3]);
            ylim([.5 .9]-(imeas-1)*[.5 .3]);
            xticks(linspace(.5-(imeas-1)*.5,.9-(imeas-1)*.3,5));
            yticks(linspace(.5-(imeas-1)*.5,.9-(imeas-1)*.3,5));
            xlabel('test');
            ylabel('retest');
            set(gca,'TickDir','out');
            set(gca,'FontSize',14);
        end
    end
end

%% Reversal Curves (accuracy / p(choosing post state)) 
trl_len = 12;
t_side = 6;

epi_start = [1 7 13 19];

p_s2     = nan(nsubj,t_side*2,2); % data
if comparesim; p_s2_sim = nan(nsubj,t_side*2,2); end % simulation

for isubj = 1:nsubj
    if ~idx_subj(isubj)
        continue
    end
    
    for icond = 0:1 % 0: bandit, 1: apples
        is_s2_temp = nan(nepis,trl_len);
        if comparesim; is_s2_sim_temp = nan(nepis*nsim,trl_len); end

        % organize data for reversal curves
        for iepi = 1:nepis
            pointer = find(idx_epi(isubj,:) == iepi & idx_cond(isubj,:) == icond,1,'first'); % index of trial after reversal 
            
            if ismember(iepi,epi_start)
                idx = pointer + (0:(t_side-icond));
            else
                idx = pointer + (-t_side:(t_side-icond)); % bandit tasks get 0 minus, apple tasks get 1
            end
            % subjects
            temp_corr = idx_blmn(isubj,idx) == idx_bmstate(isubj,idx);
            % simulations
            if comparesim
                temp_corr_sim = sim_resps(:,idx,isubj) == repmat(idx_bmstate(isubj,idx),[nsim 1]);
            end

            % if a non-starter episode, flip the first half
            if ~ismember(iepi,epi_start) 
                temp_corr(1:t_side) = 1 - temp_corr(1:t_side);
                if comparesim
                    temp_corr_sim(:,1:t_side) = 1 - temp_corr_sim(:,1:t_side);
                end
            end
            % shift bandit task for alignment w/ fairy
            if icond == 0
                temp_corr = temp_corr(2:end);
                if comparesim; temp_corr_sim = temp_corr_sim(:,2:end); end
            end
            % align initial episodes w/ the rest
            if ismember(iepi,epi_start)
                is_s2_temp(iepi,(t_side+1):end) = temp_corr;
                if comparesim; is_s2_sim_temp((1:nsim)+((iepi-1)*nsim),t_side+1:end) = temp_corr_sim; end
            else
                is_s2_temp(iepi,:) = temp_corr;
                if comparesim; is_s2_sim_temp((1:nsim)+((iepi-1)*nsim),:) = temp_corr_sim; end
            end
        end
        % calculate accuracy
        p_s2(isubj,:,icond+1) = nanmean(is_s2_temp,1);
        if comparesim; p_s2_sim(isubj,:,icond+1) = nanmean(is_s2_sim_temp,1); end;
    end
end

rgbs = [93 74 25; 25 42 68]/100*255;
% config data for plot
dat.X1 = p_s2;
dat.dlabel = {'data'};
if comparesim
    dat.X2 = p_s2_sim;
    dat.dlabel{end+1} = 'sim';
end
dat.revtrl = 7;
dat.rgbs = rgbs;
dat.glabel = condstr;
dat.figprops.YTicks = 0:.1:1;
dat.figprops.XLabel = 'trials around reversal';
dat.figprops.YLabel = 'p(post-state)';

% plot
plotReversalCurves(dat);

%% Reversal Curves (switches)
trl_len = 12;
t_side = 6;
epi_excl = [1 7 13 19];

prepe = nan(nsubj,trl_len-1,2);
if comparesim; prepe_sim = nan(nsubj,trl_len-1,2); end

for isubj = 1:nsubj
    if ~idx_subj(isubj)
        continue
    end
    
    for icond = 0:1
        repe_temp = nan(nepis-numel(epi_excl),trl_len-1);
        if comparesim; repe_temp_sim = nan((nepis-numel(epi_excl))*nsim,trl_len-1); end
        for iepi = setdiff(1:nepis,epi_excl)
            pointer = find(idx_epi(isubj,:) == iepi & idx_cond(isubj,:) == icond,1,'first'); % index of trial after reversal 

            idx = pointer + (-t_side:(t_side-icond)); % bandit tasks get 0 minus, apple tasks get 1
                
            % calculate repetitions
            temp_repe = idx_blmn(isubj,idx(2:end)) == idx_blmn(isubj,idx(1:end-1));
            if comparesim
                temp_repe_sim = sim_resps(:,idx(2:end),isubj) == sim_resps(:,idx(1:end-1),isubj);
            end

            if icond == 0
                temp_repe = temp_repe(2:end);
                if comparesim; temp_repe_sim = temp_repe_sim(:,2:end); end
            end
            repe_temp(iepi,:) = temp_repe;
            if comparesim; repe_temp_sim((1:nsim)+((iepi-1)*nsim),:) = temp_repe_sim; end
        end
        prepe(isubj,:,icond+1) = nanmean(repe_temp,1);
        if comparesim; prepe_sim(isubj,:,icond+1) = nanmean(repe_temp_sim,1); end
    end
end

rgbs = [93 74 25; 25 42 68]/100*255;

% config data for plot
dat.X1 = 1-prepe;
dat.dlabel = {'data'};
if comparesim
    dat.X2 = 1-prepe_sim;
    dat.dlabel{end+1} = 'sim';
end
dat.revtrl = 6;
dat.rgbs = rgbs;
dat.glabel = condstr;
dat.figprops.YTicks = 0:.1:1;
dat.figprops.XLabel = 'trials around reversal';
dat.figprops.YLabel = 'p(switch)';

% plot
plotReversalCurves(dat);

%% Repetition curves (signed evidence) with bm_state flags
% bandit: normal repetition curves as a function of reward
% apples: same as above but the evidence is signed 

% bin ranges
bin_edge_right = -.5:.1:.5; % previous feedback, current orange-ness

prepeat = nan(nsubj,numel(bin_edge_right)-1);
porange = nan(nsubj,numel(bin_edge_right)-1);

if comparesim
    prepeat_sim = nan(nsubj,numel(bin_edge_right)-1);
    porange_sim = nan(nsubj,numel(bin_edge_right)-1);
end

for isubj = 1:nsubj
    if ~idx_subj(isubj)
        continue
    end
    
    for icond = 0:1 % 0: bandit, 1: apples
        switch icond
            case 0 % bandit
                fb_seen = [];
                repeats = [];
                if comparesim
                    fb_seen_sim = [];
                    repeats_sim = []; 
                end

                for iblk = 1:nblk
                    idx     = idx_blk(isubj,:) == iblk & idx_cond(isubj,:) == icond; % match block and condition number
                    % data
                    fb_blk  = idx_fb(isubj,idx);
                    fb_seen = [fb_seen fb_blk(1:end-1)]; % last feedback doesn't matter
                    choices_blk = idx_blmn(isubj,idx);
                    repeats     = [repeats choices_blk(2:end) == choices_blk(1:end-1)];
                    
                    % simulations
                    if comparesim
                        fb_blk_sim  = idx_fb_sim(:,idx,isubj);
                        fb_seen_sim = [fb_seen_sim fb_blk_sim(:,1:end-1)];
                        choices_blk_sim = sim_resps(:,idx,isubj);
                        repeats_sim     = [repeats_sim choices_blk_sim(:,2:end) == choices_blk_sim(:,1:end-1)];
                    end
                end
                fb_seen = fb_seen/100-.5; % data
                if comparesim; fb_seen_sim = fb_seen_sim/100-.5; end % sim
                
                for ibin = 2:numel(bin_edge_right)
                    idx_binrange = fb_seen > bin_edge_right(ibin-1) & fb_seen <= bin_edge_right(ibin);
                    prepeat(isubj,ibin-1) = sum(repeats(idx_binrange))/sum(idx_binrange);
                    
                    if comparesim
                        idx_binrange = fb_seen_sim > bin_edge_right(ibin-1) & fb_seen_sim <= bin_edge_right(ibin);
                        prepeat_sim(isubj,ibin-1) = sum(repeats_sim(idx_binrange))/sum(idx_binrange,'all');
                    end
                end
                
            case 1 % apples
                fb_signed   = [];
                repeats     = [];
                if comparesim
                    fb_signed_sim = [];
                    repeats_sim = []; 
                end

                for iblk = 1:nblk
                    idx     = idx_blk(isubj,:) == iblk & idx_cond(isubj,:) == icond;
                    % data
                        fb_blk  = idx_fb(isubj,idx)/100-.5; % fb_blk is the strength of the evidence FOR the current state
                        % make blue the positively signed absolute evidence
                        idx_blue = idx_bmstate(isubj,idx) == 1;
                        fb_blk(~idx_blue) = 1-fb_blk(~idx_blue);
                        fb_blk  = fb_blk(2:end);% the 1st trial of a block has nothing to be signed with
                        % calculate sign vector
                        signer = idx_blmn(isubj,idx);
                        signer(signer == 0) = -1;
                        % sign evidence at t with choice from t-1
                        fb_blk = fb_blk .* signer(1:end-1);
                        fb_signed = [fb_signed fb_blk];
                        % calculate repeats
                        choices_blk = idx_blmn(isubj,idx);
                        repeats = [repeats choices_blk(2:end) == choices_blk(1:end-1)];

                    % sim
                    if comparesim
                        fb_blk_sim = idx_fb_sim(:,idx,isubj)/100-.5;
                        % make blue the positively signed absolute evidence
                        idx_blue = repmat(idx_bmstate(isubj,idx),[nsim 1]) == 1;
                        fb_blk_sim(~idx_blue) = 1-fb_blk_sim(~idx_blue);
                        fb_blk_sim(:,1)  = []; % the 1st trial of a block has nothing to be signed with
                        % calculate sign vector
                        signer = sim_resps(:,idx,isubj);
                        signer(signer == 0) = -1;
                        % sign evidence at t with choice from t-1
                        fb_blk_sim = fb_blk_sim .* signer(:,1:end-1);
                        fb_signed_sim = [fb_signed_sim fb_blk_sim];
                        % calculate repeats
                        choices_blk_sim = sim_resps(:,idx,isubj);
                        repeats_sim = [repeats_sim choices_blk_sim(:,2:end) == choices_blk_sim(:,1:end-1)];
                    end
                    
                end
                
                for ibin = 2:numel(bin_edge_right)
                    idx_binrange = fb_signed > bin_edge_right(ibin-1) & fb_signed <= bin_edge_right(ibin);
                    porange(isubj,ibin-1) = sum(repeats(idx_binrange))/sum(idx_binrange);
                    if comparesim
                        idx_binrange = fb_signed_sim > bin_edge_right(ibin-1) & fb_signed_sim <= bin_edge_right(ibin);
                        porange_sim(isubj,ibin-1) = sum(repeats_sim(idx_binrange))/sum(idx_binrange,'all');
                    end
                end
        end
    end
end

xlen = numel(bin_edge_right)-1;
ndat = sum(any(~isnan(porange(:,:,1)),2));
rgbs = [93 74 25; 25 42 68]/100*255;

% config data for plot
dat.X1 = cat(3,prepeat,porange);
dat.dlabel = {'data'};
if comparesim
    dat.X2 = cat(3,prepeat_sim,porange_sim);
    dat.dlabel{end+1} = 'sim';
end
dat.rgbs = rgbs;
dat.glabel = condstr;
dat.figprops.YTicks = 0:.1:1;
dat.figprops.XLabel = 'bins (signed feedback)';
dat.figprops.YLabel = 'p(repeat)';

% plot
plotReversalCurves(dat);
plot([1 xlen],[0 1],'--','HandleVisibility','off');
title(sprintf('p(repeat) - Bandit\np(repeat, ev. signed by prev choice) - Inference\n%s (N = %d)',samplename,nsubj-nexcl),'FontSize',14);
yline(.5,':','HandleVisibility','off');
xline(5.5,':','HandleVisibility','off');
xticklabels((bin_edge_right(2:end)+bin_edge_right(1:end-1))/2);
xlabel('bins (signed feedback)');
hold off


%% Logistic regression

nreg = 5; % number of regressors for logistic regression analysis

breg_cond = nan(nsubj,nreg+1,2);
pval_cond = nan(nsubj,nreg+1,2);
for isubj = 1:nsubj
    if ~idx_subj(isubj)
        continue
    end
    for icond = 1:2
        % get current subject and condition
        subj = isubj;
        cond = icond-1;

        % get task variables and behavior
        ifilt = idx_cond(subj,:) == cond;

        blknum  = idx_blk(subj,ifilt);
        trlnum  = repmat(1:nnz(blknum == 1),[1,numel(unique(blknum))]);
        bmstate = 2-idx_bmstate(subj,ifilt);
        fbcorr  = idx_fbabs(subj,ifilt)/100*2-1;
        fbboth  = cat(1,+fbcorr,-fbcorr);
        xevi    = fbboth(sub2ind(size(fbboth),bmstate,1:numel(bmstate)))*xfac;
        resp    = 2-idx_blmn(subj,ifilt);
        iscor   = idx_corr(subj,ifilt);
        clear fbcorr fbboth

        % realign the two conditions on a common frame
        % response (resp) now follows the evidence (xevi) in both conditions
        if cond == 0 % bandit
            itrl = find(trlnum <= 72);
            blknum = blknum(itrl);
            trlnum = trlnum(itrl);
            bmstate = bmstate(itrl);
            xevi = xevi(itrl);
            resp = resp(itrl+1);
        else % fairy
            itrl = find(trlnum >= 2);
            blknum = blknum(itrl);
            trlnum = trlnum(itrl)-1;
            bmstate = bmstate(itrl);
            xevi = xevi(itrl);
            resp = resp(itrl);
        end

        % run logistic regression analysis
        itrl = find(trlnum > nreg);
        xreg = xevi(itrl)';
        for ireg = 1:nreg
            xreg = cat(2,xreg,xevi(itrl-ireg)');
        end
        yreg = resp(itrl)' == 1;
        [breg,~,stat] = glmfit(xreg,yreg,'binomial','link','logit');
        breg_cond(isubj,:,icond) = breg(2:end);
        pval_cond(isubj,:,icond) = stat.p(2:end);
    end
end

% fraction of subjects not using evidence wrt its position
p0 = squeeze(sum(pval_cond > 0.05,1));
figure;
hold on
xlim([-0.5,nreg+0.5]);
ylim([0,1]);
plot(xlim,[0.5,0.5],'k--');
plot(0:nreg,p0/nsubj,'-','LineWidth',2);
hold off
xlabel('evidence position from current trial');
ylabel('fraction of subjects not using evidence');

% regression coefficient of evidence wrt its position
bq = squeeze(quantile(breg_cond,[0.25,0.50,0.75],1));
figure;
for icond = 1:2
    subplot(1,2,icond);
    hold on
    ylim([-0.5,3]);
    bar(0:nreg,bq(2,:,icond));
    for i = 1:size(bq,2)
        plot([i,i]-1,bq([1,3],i,icond),'k-');
    end
    hold off
    xlabel('evidence position from current trial');
    ylabel('regression coefficient');
end


%% Retest similarity
rt_simil = nan(nsubj,ntrl*4,ncnd);

for isubj = 1:nsubj
    if ~idx_rt(isubj)
        continue
    end
    for icond = 0:1
        for iblk = 1:4
            idx_t  = idx_cond(isubj,:) == icond & idx_blk(isubj,:) == iblk;
            idx_r = dat_rt.idx_cond(isubj,:) == icond & dat_rt.idx_blk(isubj,:) == iblk;
            
            rt_simil(isubj,(1+73*(iblk-1)):(73*(iblk)),icond+1) = dat_rt.idx_blmn(isubj,idx_r) == idx_blmn(isubj,idx_t);
        end
    end
end

