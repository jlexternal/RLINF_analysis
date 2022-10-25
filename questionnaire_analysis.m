% modelfree_analysis
clear all

addpath ./plot_functions/
samplename = 'pilot07'; % Pilot - 'pilot'
load(sprintf('./processed/%s/preprocessed_data_%s.mat',samplename,samplename)); % load the raw data structure for data sample
load(sprintf('./processed/%s/ques_struct.mat',samplename));
ques_excl = load(sprintf('./processed/%s/idx_excl_ques.mat',samplename));

nsubj = size(idx_fb,1);
idx_excl_task = isnan(idx_fb(:,1));
nexcl = sum(idx_excl_task);
idx_excl_ques = false(nsubj,1);
idx_excl_ques(ques_excl.idx_excl_ques) = 1;

idx_excl_all = idx_excl_task | idx_excl_ques;

fprintf('Number of subjects: %d\n',nsubj-nexcl);
condstr = {'bandit','apples'};

ques_order = {'depress', 'anxiety', 'ocir', 'social', 'bis', 'schizo', 'alcohol', 'eat', 'apathy', 'iq'};

% load questionnaire scores into compact structure
for isubj = 1:nsubj
    if idx_excl_all(isubj)
        continue
    end
    ques = ques_struct{isubj};

    % initialize questionnaire score structure
    if ~exist('q_scores','var')
        q_scores = struct;
        q_names = fields(ques)';
        for name = q_names
            q_scores.(name{1}) = nan(1,nsubj);
        end
    end

    % load in subject questionnaire scores into structure
    for q = fields(ques)'
        qname = q{1};
        if isstruct(ques.(qname).score)
            q_scores.(qname)(isubj) = ques.(qname).score.total;
        else
            q_scores.(qname)(isubj) = ques.(qname).score;
        end
        % error check (nan assignment)
        if isnan(q_scores.(qname)(isubj))
            error('error found with subj %d! (NaN assignment on included subject) \n',isubj)
        end
    end
end

q_scores = orderfields(q_scores, ques_order);

%% histogram of questionnaire scores

load ./constants/constants_ques.mat % loads ques_order, min_ques, max_ques, quesrgb

figure(1);
clf
for iques = 1:numel(ques_order)
    subplot(2,5,iques)
    histogram(q_scores.(ques_order{iques}));
    xline(min_ques(iques));
    xline(max_ques(iques));

end





