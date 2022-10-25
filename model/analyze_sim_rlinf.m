% analyze_sim_rlinf
clear all

addpath ../../toolbox/plot_functions/

samplename = 'pilot07';

% load data
load(sprintf('../../processed/%s/preprocessed_data_%s.mat',samplename,samplename));
load(sprintf('../../constants/constants_rlinf_%s.mat',samplename)); % load nblk, ncond, ntrl, samplename

nsubj = size(idx_blmn,1);
nexcl = sum(isnan(idx_blmn(:,1)));
npar = 3;
ncond = 2;
condstr = {'bandit','apples'};

idx_subj = true(nsubj,1) & ~isnan(idx_blmn(:,1));

% --------- (SPECIFY FILENAME HERE) -----------
load ./out/out_sim_nsim1000_pilot07.mat % load simulation data

if ~strcmpi(sim_out.samplename,samplename)
    error('Sample ID of data does not match that of simulations!');
end
% matrix holding simulation responses
nsim = 1e3; %sim_out.nsim;
resp_sim = nan(nsim,ntrl*nblk*ncond,nsubj);

%% Reversal curves: Accuracy / P(post-state)
t_len = 12;
t_side = 6;
p_s2 = nan(nsubj,t_len,ncond);

epi_start = [1 7 13 19];

for isubj = 1:nsubj
    if isnan(idx_blmn(isubj,1))
        continue
    end
    % extract simulation responses in the same form as subjects
    resp_sim(:,:,isubj) = cat(2,sim_out.out{isubj,1}.rt,sim_out.out{isubj,2}.rt);
    is_bmstate = idx_bmstate(isubj,:);
    is_bmstate(is_bmstate == 0) = 2; % convert isbmstate to 1/2 from 1/0
    corr_sim_temp = resp_sim(:,:,isubj) == is_bmstate;

    for icond = 0:1
        is_s2_temp = nan(nepis*nsim,t_len); % flag for post-state response
        for iepi = 1:nepis
            % pointer at trial post-reversal
            ptr = find(idx_epi(isubj,:) == iepi & idx_cond(isubj,:) == icond,1,'first'); 

            % if the beginning of a block
            if ismember(iepi,epi_start)
                idx = ptr + (0:(t_side-icond));
            else
                idx = ptr + (-t_side:(t_side-icond));
            end
            temp_corr = corr_sim_temp(:,idx);

            % if a non-starter episode, don't flip the first half
            if ~ismember(iepi,epi_start) 
                temp_corr(:,1:t_side) = 1-temp_corr(:,1:t_side);
            end
            
            % shifting for bandit tasks compared to fairy task
            if icond == 0
                temp_corr = temp_corr(:,2:end);
            end

            if ismember(iepi,epi_start)
                is_s2_temp((1:nsim)+((iepi-1)*nsim),t_side+1:end) = temp_corr;
            else
                is_s2_temp((1:nsim)+((iepi-1)*nsim),:) = temp_corr;
            end
        end
        p_s2(isubj,:,icond+1) = nanmean(is_s2_temp,1);
    end

end

xlen = size(p_s2,2);
ndat = sum(any(~isnan(p_s2(:,:,1)),2));
rgbs = [93 74 25; 25 42 68]/100;

figure(2);
clf
for icond = 1:2
    curve = mean(p_s2(idx_subj,:,icond),1);
    curvesem = std(p_s2(idx_subj,:,icond),[],1)/sqrt(ndat);
    
    shadedErrorBar(1:xlen,curve,curvesem,'lineprops',{'Color',rgbs(icond,:),'LineWidth',2});
end
title(sprintf('Reversal curves\np(post-state)\nfrom SIMULATIONS\n%s (N = %d, nsim = %d)',samplename,nsubj-nexcl,nsim),'FontSize',14);
ylim([0 1]);
yline(.5,':','HandleVisibility','off');
xticks(1:12);
xticklabels({'-6','-5','-4','-3','-2','-1','1','2','3','4','5','6'});
xlabel('trials around reversal');
legend({'bandit','apples'},'Location','southeast','FontSize',14);
xline(6.5,'HandleVisibility','off');



