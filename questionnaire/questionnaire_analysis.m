% modelfree_analysis
clear all

addpath ../toolbox/
addpath ../toolbox/plot_functions/
% ----------------- Input -----------------
samplename = 'sample2'; % Pilot - 'pilot'
isretest = false;
% -----------------------------------------
load(sprintf('../constants/constants_rlinf_%s.mat',samplename));
rt_str = '';
if isretest
    rt_str = 'retest_';
    fprintf('Analyzing retest data...\n');
end
load(sprintf('../processed/%s/preprocessed_data_%s%s.mat',samplename,rt_str,samplename)); % load the raw data structure for data sample
load(sprintf('../processed/%s/ques_struct_%s%s.mat',samplename,rt_str,samplename));
if isretest
    load(sprintf('../processed/%s/idx_TaskQuesAll_rt.mat',samplename),'idx_fullAll_rt');
    idx_fullAll = idx_fullAll_rt;
else
    load(sprintf('../processed/%s/idx_TaskQuesAll.mat',samplename),'idx_fullAll');
end

nsubj = size(idx_blmn,1);
fprintf('Number of subjects: %d\n',sum(idx_fullAll));

ques_order = {'depress', 'anxiety', 'ocir', 'social', 'bis', 'schizo', 'alcohol', 'eat', 'apathy', 'iq'};
ques_excl = {'schizo','alcohol'};

% load raw questionnaire responses into compact structure
for isubj = 1:nsubj
    if ~idx_fullAll(isubj)
        continue
    end
    ques = ques_struct{isubj};
    % initialize questionnaire score structure
    if ~exist('q_scores','var')
        q_scores = struct;
        q_names = fields(ques)';
        for name = q_names
            nitem = numel(ques_struct{isubj}.(name{1}).raw);
            if strcmpi(name,'iq')
                nitem = 1;
            elseif strcmpi(name,'social')
                nitem = numel(ques_struct{isubj}.(name{1}).avg);
            end
            q_scores.(name{1}) = nan(nsubj,nitem);
        end
        % add excluded questionnaire for placekeeping
        q_scores.schizo  = nan(nsubj,43);
        q_scores.alcohol = nan(nsubj,10);
    end
    % load in subject questionnaire scores into structure
    for q = fields(ques)'
        qname = q{1};
        if strcmpi(qname,'iq')
            q_scores.(qname)(isubj) = ques_struct{isubj}.(qname).score;
            continue
        elseif strcmpi(qname,'social')
            q_scores.(qname)(isubj,:) = ques_struct{isubj}.(qname).avg;
            continue
        end
        q_scores.(qname)(isubj,:) = ques.(qname).raw;
    end
end
q_scores = orderfields(q_scores, ques_order);

%% histogram of questionnaire scores

load ../constants/constants_ques.mat % loads ques_order, min_ques, max_ques, quesrgb

figure(1);
clf
sgtitle(sprintf('Questionnaire itemsum distributions\nnsubj = %d',nsubj));
for iques = 1:numel(ques_order)
    subplot(2,5,iques);
    hold on
    title(sprintf('%s',ques_order{iques}));
    scores_summed = nansum(q_scores.(ques_order{iques}),2);
    if any(strcmpi(ques_order{iques},ques_excl))
        continue
    end
    histogram(scores_summed);
    xline(min(scores_summed));
    xline(max(scores_summed));
end

%% convert raw questionnaire scores into scaled z-scored values and project onto trandiagnostic scores

% 1/ Get z-score distribution of log-transformed scores of Rouault et al (2018)
load ./external/zparams_marion.mat % loads zparams
% 2/ Log-transform raw scores of RLINF
scores_zlogt = [];

for q = ques_order
    field = q{1};
    if strcmpi(field,'iq')
        continue
    elseif any(strcmpi(field,{'ocir','social','schizo','eat','alcohol'}))
        sc = 1+q_scores.(field); % to avoid blow-up of log(0)
    else
        sc = q_scores.(field);
    end
    sc_logt = log(sc);
    if any(isinf(sc_logt))
        error('Inf found in log-transform!')
    end
    % choose appropriate z-transform
    if strcmpi(field,'depress')
        fprintf('1: %s\n',field)
        zmu  = zparams.('zung').mu;
        zsig = zparams.('zung').sigma;
    elseif strcmpi(field,'social')
        fprintf('2: %s\n',field)
        zmu  = zparams.('leb').mu;
        zsig = zparams.('leb').sigma;
    else
        fprintf('3: %s\n',field)
        zmu  = zparams.(field).mu;
        zsig = zparams.(field).sigma;
    end
    
    % apply transform
    sc_zlogt = (sc_logt - zmu)/zsig; 
    if any(isinf(sc_zlogt))
        error('Inf found in z-scoring!')
    end
    scores_zlogt = cat(1,scores_zlogt,sc_zlogt');
end

isubj = find(idx_fullAll,1,'first');
idx_nan = isnan(scores_zlogt(:,isubj));
scores_zlogt(idx_nan,idx_fullAll) = 0;
% 4/ Matrix multiply mini-questionnaire coefficients to obtain t-diag scores
load ./external/coef_mini_q.mat

% calculate scores
sc_dim = nan(nsubj,3);
for idim = 1:3
    sc_dim(:,idim) = sum(coef_mini(:,idim).*scores_zlogt,1);
end

save(sprintf('../processed/%s/dim_scores_mini_%s%s.mat',samplename,rt_str,samplename),'sc_dim');
% plot
addpath ../toolbox/plot_functions/
props = struct;
props.MomentType = 'median';
props.ErrorType = 'quartiles';
plotViolins(sc_dim(idx_fullAll,:),psychrgb*255,props,psychstr);

%% test-retest reliability of mini t-diag scores

% load appropriate variables here (to-be-finalized)
sc_dim = load(sprintf('../processed/%s/dim_scores_mini_%s.mat',samplename,samplename));
sc_dim = sc_dim.sc_dim;
sc_dim_rt = load(sprintf('../processed/%s/dim_scores_mini_retest_%s.mat',samplename,samplename));
sc_dim_rt = sc_dim_rt.sc_dim;
for idim = 1:3
    x = sc_dim(idx_fullAll_rt,idim);
    y = sc_dim_rt(idx_fullAll_rt,idim);
    xrange = [min(x) max(x)];
    [pn,s] = polyfit(x,y,1);
    [py,d] = polyconf(pn,xrange,s,'alpha',0.05,'predopt','curve');
    figure(idim);
    clf
    hold on
    s = shadedErrorBar(xrange,py,d,'patchSaturation',.1,'lineprops',{'LineWidth',1,'Color',psychrgb(idim,:)});
    set(s.edge,'LineStyle','none');
    scatter(x,y,'MarkerFaceColor',psychrgb(idim,:),'MarkerEdgeColor','none','MarkerFaceAlpha',.5);
    [r,p] = corr(x,y);
    plot([-2 3],[-2 3],':');
    xlim([min(x)-.2 max(x)+.2]);
    ylim([min(y)-.2 max(y)+.2]);
    xlabel(sprintf('%s (test)',psychstr{idim}));
    ylabel(sprintf('%s (retest)',psychstr{idim}))
    title(sprintf('%s (test-retest)\nr^2=%.4f',psychstr{idim},r^2),'FontSize',12);

    fprintf('Test-retest ICC(%s):\n',psychstr{idim});
    fprintf('p(correct): Pearson r= %.4f (p=%.4f)\n',r,p);
    [icc_hat,icc_loc,icc_hic,~,~,~,pval] = ICC([x,y],'C-1');
    fprintf('p(correct): ICC(3,1) = %.3f [%.3f %.3f] (p=%.4f)\n',icc_hat,icc_loc,icc_hic,pval);
    [icc_hat,icc_loc,icc_hic,~,~,~,pval] = ICC([x,y],'1-1');
    fprintf('p(correct): ICC(1,1) = %.3f [%.3f %.3f] (p=%.4f)\n',icc_hat,icc_loc,icc_hic,pval);
    end



