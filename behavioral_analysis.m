% behavioral_analysis
clear all

addpath ./toolbox/plot_functions/
samplename = 'sample2'; % Pilot - 'pilot'

load(sprintf('./constants/constants_rlinf_%s.mat',samplename)); % load nblk, ncond, ntrl, samplename
load(sprintf('./processed/%s/preprocessed_data_%s.mat',samplename,samplename)); % load the raw data structure for data sample

nsubj = size(idx_fb,1);
nexcl = sum(isnan(idx_fb(:,1)));

fprintf('Number of subjects: %d\n',nsubj-nexcl);
condstr = {'bandit','fairy'};

idx_firsttrl = 1+73*(0:7);

% -------------------------INPUT--------------------------------
% Compare participant behavior with simulations
comparesim = false; 
if comparesim
    simfilename = 'processed_out_sim_noisyINF_nsim1000_pilot07.mat';
    load(sprintf('./processed/%s/%s',samplename,simfilename));
    sim_resps = out.resps;
    sim_resps(sim_resps == 2) = 0;
    nsim = out.nsim;
end
% --------------------------------------------------------------

%% Overall accuracy and repetitions
pcorr = nan(nsubj,2);
prepe = nan(nsubj,2);

for isubj = 1:nsubj
    if isnan(idx_fb(isubj,1))
        continue
    end
    
    for icond = 0:1
        % Calculating accuracy: ignore the 1st trial of each block
        idx = idx_cond(isubj,:) == icond & ~ismember(idx_abstr(isubj,:),idx_firsttrl);
        pcorr(isubj,icond+1)    = nanmean(idx_corr(isubj,idx));
        
        % repetitions must be calculated within each block
        nrepe = 0;
        for iblk = 1:4
            idx = idx_cond(isubj,:) == icond & idx_blk(isubj,:) == iblk;
            resps = idx_blmn(isubj,idx);
            if icond == 0 % bandit (compare starting 2nd to 1st trial)
                nrepe = nrepe + sum(resps(2:end) == resps(1:end-1));
            else % apples: (compare starting 3rd to 2nd trial, since first trial is nothing)
                nrepe = nrepe + sum(resps(3:end) == resps(2:end-1));
            end
        end
        
        if icond == 0
            ntot = 72;
        else
            ntot = 71;
        end 
        prepe(isubj,icond+1) = nrepe/(ntot*4);
    end
end

if true
    save(sprintf('./processed/%s/pcor_%s.mat',samplename,samplename),'pcorr');
end

fprintf('General performance (%s, N=%d)\n',samplename,nsubj-nexcl);
for icond = 1:2
    fprintf('%s\n',condstr{icond});
    fprintf('Accuracy: %.2f ± %.2f\n',nanmean(pcorr(:,icond)),nanstd(pcorr(:,icond)));
    fprintf('Repeats: %.2f ± %.2f\n\n',nanmean(prepe(:,icond)),nanstd(prepe(:,icond)));
end

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

%% Reversal Curves (accuracy / p(choosing post state)) 
trl_len = 12;
t_side = 6;

epi_start = [1 7 13 19];

p_s2     = nan(nsubj,t_side*2,2); % data
if comparesim; p_s2_sim = nan(nsubj,t_side*2,2); end % simulation

for isubj = 1:nsubj
    if isnan(idx_fb(isubj,1))
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
    if isnan(idx_fb(isubj,1))
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

%% repetition curves (signed evidence) with bm_state flags
% bandit: normal repetition curves as a function of reward
% apples: same as above but the evidence is signed 

% bin ranges
bin_edge_right = -.5:.1:.5; % previous feedback, current orange-ness

prepeat = nan(nsubj,numel(bin_edge_right)-1);
porange = nan(nsubj,numel(bin_edge_right)-1);

for isubj = 1:nsubj
    if isnan(idx_fb(isubj,1))
        continue
    end
    
    for icond = 0:1 % 0: bandit, 1: apples
        switch icond
            case 0 % bandit
                fb_seen = [];
                repeats = [];
                for iblk = 1:nblk
                    idx     = idx_blk(isubj,:) == iblk & idx_cond(isubj,:) == icond; % match block and condition number
                    fb_blk  = idx_fb(isubj,idx);
                    fb_seen = [fb_seen fb_blk(1:end-1)]; % last feedback doesn't matter
                    
                    choices_blk = idx_blmn(isubj,idx);
                    repeats     = [repeats choices_blk(2:end) == choices_blk(1:end-1)];
                end
                
                fb_seen = fb_seen/100-.5;
                for ibin = 2:numel(bin_edge_right)
                    idx_binrange = fb_seen > bin_edge_right(ibin-1) & fb_seen <= bin_edge_right(ibin);
                    prepeat(isubj,ibin-1) = sum(repeats(idx_binrange))/sum(idx_binrange);
                end
                
            case 1 % apples
                fb_signed   = [];
                repeats     = [];
                for iblk = 1:nblk
                    idx     = idx_blk(isubj,:) == iblk & idx_cond(isubj,:) == icond;
                    
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
                end
                
                for ibin = 2:numel(bin_edge_right)
                    idx_binrange = fb_signed > bin_edge_right(ibin-1) & fb_signed <= bin_edge_right(ibin);
                    porange(isubj,ibin-1) = sum(repeats(idx_binrange))/sum(idx_binrange);
                end
        end
    end
end

xlen = numel(bin_edge_right)-1;
ndat = sum(any(~isnan(porange(:,:,1)),2));
rgbs = [93 74 25; 25 42 68]/100;
figure
clf
hold on
shadedErrorBar(1:xlen,nanmean(prepeat,1),nanstd(prepeat,[],1)/sqrt(ndat),'lineprops',{'Color',rgbs(1,:),'LineWidth',2});
shadedErrorBar(1:xlen,nanmean(porange,1),nanstd(porange,[],1)/sqrt(ndat),'lineprops',{'Color',rgbs(2,:),'LineWidth',2});
plot([1 xlen],[0 1],'--');
title(sprintf('p(repeat) - Bandit\np(repeat, ev. signed by prev choice) - Inference\n%s (N = %d)',samplename,nsubj-nexcl),'FontSize',14);
xticklabels((bin_edge_right(2:end)+bin_edge_right(1:end-1))/2);
xlabel('bins (signed feedback)');
legend(condstr,'Location','southeast','FontSize',14);
yline(.5,':','HandleVisibility','off');
xline(5.5,':','HandleVisibility','off');
ylim([0 1]);
xlim([.5 10.5]);

%% choice similarity matrices 

emax = 25;
eavg = 12;
emin = 6;

nepi = 24;

ntotal = zeros(emax,emax,nsubj,2); % last index is condition
nrepeat = zeros(emax,emax,nsubj,2);

for isubj = 1:nsubj
    if isnan(idx_fb(isubj,1))
        continue
    end
    
    for icond = 0:1
        for iepi = 1:nepi
            idx = idx_cond(isubj,:) == icond & idx_epi(isubj,:) == iepi;
            nt = sum(idx);
            
            % convert responses to complex vector
            z_resp = double(idx_blmn(isubj,idx) == idx_blmn(isubj,find(idx==1,1,'first')));
            z_resp = complex(z_resp, double(~z_resp));
            
            nrepeat(1:nt,1:nt,isubj,icond+1) = nrepeat(1:nt,1:nt,isubj,icond+1) + real(z_resp'.*z_resp);
            ntotal(1:nt,1:nt,isubj,icond+1)  = ntotal(1:nt,1:nt,isubj,icond+1) + ones(nt,nt);
        end
    end
end

nt_plot = eavg+4;

psim = nrepeat(1:nt_plot,1:nt_plot,:,:)./ntotal(1:nt_plot,1:nt_plot,:,:);

figure(4)
clf
for icond = 1:2
    subplot(1,2,icond);
    title(sprintf('%s',condstr{icond}),'FontSize',14);
    hold on
    colormap('hot');
    imagesc(triu(nanmean(psim(:,:,:,icond),3),0));
    if icond == 2
        colorbar;
    end
    
    xlim([0,nt_plot]+.5);
    ylim([0,nt_plot]+.5);
    xticks(0:4:nt_plot);
    yticks(0:4:nt_plot);
    set(gca,'TickDir','out','FontSize',14);
end
sgtitle(sprintf('choice similarity matrices\n%s (N = %d)',samplename,nsubj-nexcl),'FontSize',14);

