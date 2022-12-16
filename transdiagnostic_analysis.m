% transdiagnostic_analysis
clear all
clc
% ------------------------ INPUT: ------------------------
samplename = 'sample2';
modelkernel = 'noisyKF';
% --------------------------------------------------------
addpath(genpath('./toolbox/cprintf'));
% load constants
load(sprintf('./constants/constants_rlinf_%s.mat',samplename));
% load parameter fits from chosen model
load(sprintf('./frontex_environment/import/sample_out/%s/pars_fit_%s_ALL.mat',samplename,modelkernel),'out');
pars = out.pars;
xnam = out.xnam;
fprintf('Loading parameters from model with description:\n ''%s''\n',out.description);
clearvars out
try
    load(sprintf('./frontex_environment/import/sample_out/%s/pars_fit_%s_retest_ALL.mat',samplename,modelkernel),'out');
    pars_rt = out.pars;
    clearvars out
catch
    warning('Parameter fits for %s retest cannot be found!',modelkernel);
end

% load indices 
load(sprintf('./processed/%s/idx_TaskQuesAll.mat',samplename),'idx_fullAll');
load(sprintf('./processed/%s/idx_TaskQuesAll_rt.mat',samplename),'idx_fullAll_rt');
load(sprintf('./processed/%s/idx_binomialTestPass.mat',samplename));
idx_binomialTestPass = any(idx_binomialTestPass,2);

% load transdiagnostic dimension scores
load(sprintf('./processed/%s/dim_scores_mini_%s.mat',samplename,samplename),'sc_dim');
sc_dim_rt = load(sprintf('./processed/%s/dim_scores_mini_retest_%s.mat',samplename,samplename),'sc_dim');
sc_dim(:,:,2) = sc_dim_rt.sc_dim;
clearvars sc_dim_rt

% load ICAR 
load(sprintf('./processed/%s/icar_sample2.mat',samplename),'icar');
icar_rt = load(sprintf('./processed/%s/icar_retest_sample2.mat',samplename),'icar');
icar(:,2) = icar_rt.icar;
clearvars icar_rt

% load performance measures
load(sprintf('./processed/%s/pcor_%s.mat',samplename,samplename),'pcorr');
load(sprintf('./processed/%s/pcor_rt_%s.mat',samplename,samplename),'pcorr_rt');
load(sprintf('./processed/%s/prep_%s.mat',samplename,samplename),'prepe');
load(sprintf('./processed/%s/prep_rt_%s.mat',samplename,samplename),'prepe_rt');

% choose subset of participants
%     idx_fullAll: Participants in test who completed both task and questionnaire
%     idx_binomialTestPass: Participants who passed binomial test for accuracy
%     idx_fullAll_rt: Participants in retest who completed both task and questionnaire
idx_subj = idx_fullAll & idx_binomialTestPass;

teststr = {'test','retest'};

%% Correlations of measures with transdiagnostic dimensions
corrtype = 'Pearson';

clc
fprintf('Model kernel: %s\n',modelkernel);
fprintf('Correlation type: %s\n',corrtype);
for itest = 1
    fprintf('Dataset: %s (N = %d)\n',teststr{itest},sum(idx_subj));
    if itest == 1
        % build matrix of measures
        y = [pcorr 1-prepe squeeze(pars(:,:,1)) squeeze(pars(:,:,2))];
    else % retest
        y = [pcorr_rt prepe_rt squeeze(pars_rt(:,:,1)) squeeze(pars_rt(:,:,2))];
    end
    y_str = {'p(corr)_B','p(corr)_F','p(swi)_B','p(swi)_F'};
    y_str = [y_str strcat(xnam,'_B') strcat(xnam,'_F')];

    if size(y,2) ~= numel(y_str)
        error('missing labels!')
    end
    y = y(idx_subj,:);
    for idim = 1:3
        fprintf('%s\n',psychstr{idim});
        x = sc_dim(idx_subj,idim,itest);
        [r,p] = corr(x,y,'Type',corrtype);
        for iy = 1:numel(y_str)
            if p(iy) < 0.05
                cprintf('*white',' ~ %s: r=%+.4f (p=%.4f)\n',pad(y_str{iy},10),r(iy),p(iy));
            else
                fprintf(' ~ %s: r=%+.4f (p=%.4f)\n',pad(y_str{iy},10),r(iy),p(iy));
            end
        end
    end
end

%% Regressions of measures with transdiagnostic dimensions

fprintf('Model kernel: %s\n',modelkernel);
for itest = 1
    fprintf('Dataset: %s (N = %d)\n',teststr{itest},sum(idx_subj));
    if itest == 1
        % build matrix of measures
        y = [pcorr 1-prepe squeeze(pars(:,:,1)) squeeze(pars(:,:,2))];
    else % retest
        y = [pcorr_rt prepe_rt squeeze(pars_rt(:,:,1)) squeeze(pars_rt(:,:,2))];
    end
    y_str = {'p(corr)_B','p(corr)_F','p(swi)_B','p(swi)_F'};
    y_str = [y_str strcat(xnam,'_B') strcat(xnam,'_F')];
    if size(y,2) ~= numel(y_str)
        error('missing labels!')
    end
    y = y(idx_subj,:);

    x = sc_dim(idx_subj,:,itest);
    for imeas = 1:numel(y_str)
        fprintf('%s: \n',y_str{imeas});
        mdl = fitglm(x,y(:,imeas),'VarNames',[psychstr y_str{imeas}])
        disp(' ')
    end
end




