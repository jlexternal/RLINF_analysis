% modelfree_analysis
clear all

addpath ../toolbox/plot_functions/
samplename = 'sample2'; % Pilot - 'pilot'
load(sprintf('../processed/%s/preprocessed_data_%s.mat',samplename,samplename)); % load the raw data structure for data sample
load(sprintf('../processed/%s/ques_struct.mat',samplename));
load(sprintf('../processed/%s/idx_TaskQuesAll.mat',samplename),'idx_fullAll');

nsubj = size(idx_blmn,1);
fprintf('Number of subjects: %d\n',sum(idx_fullAll));
condstr = {'bandit','apples'};

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

% plot
trgb  = [182 137 115; 234 191 159; 137 137 137];
tstr = {'AD','CIT','SW'};

addpath ../toolbox/plot_functions/
props = struct;
props.MomentType = 'mean';
props.ErrorType = 'sem';
plotViolins(sc_dim,trgb,props,{'AD','CIT','SW'});





