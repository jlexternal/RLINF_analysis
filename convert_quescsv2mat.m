% convert_quescsv2mat.m
%
% Usage: Converts questionnaire csv files taken from the MySQL database 
% (with predefined selection clauses) into MATLAB structures.
clear all
clc

addpath ./toolbox/ques_functions/ % for parseqn.m

% --------------------INPUTS ------------------------------
samplename  = 'sample2';
versionadd  = '';
nsubj       = 247; % can be arbitrary but set to constant later
dirname     = './import';
filename    = sprintf('ques_%s%s.csv',samplename,versionadd);
is_short_ques = true; % toggle for shortened questionnaire
% ---------------------------------------------------------

if is_short_ques
    nques = 11;
    warning('Organizing questionnaire data based on shortened optimized indexing!');
else
    nques = 14;
    warning('Organizing questionnaire data based on full battery indexing!');
end

test = parseqn(sprintf('./%s/%s/%s',dirname,samplename,filename),inf); % when the table ORDER BY is by unique_id and label (taken from ques_table_copy)
test(1,:) = [];

if is_short_ques
    labelstr = {'alcohol','anxiety','apathy','bis','depress','eat','icar','ocir','schizo','social'};
    labelcnt = [0 1 1 1 1 1 3 1 0 2]; % number of rows per questionnaire label
else
    labelstr = {'alcohol','anxiety','apathy','bis','depress','eat','icar','ocir','schizo','social'};
    labelcnt = [1 1 1 1 1 1 3 1 2 2]; % number of rows per questionnaire label
end

% index of included questionnaires for the shortened optimized questionnaire battery
ques_opt_incl = struct;
ques_opt_incl.bis       = [6,9,13,14,17,19,20,22,25];
ques_opt_incl.ocir      = [2,3,4,6,7,9,11,12,13,15,16,18];
ques_opt_incl.schizo    = [];
ques_opt_incl.depress   = [10,11,12,13,14,15,17,18,20];
ques_opt_incl.social    = [2,5,6,7,8,9,10,11,12,14,15,16,18,19,20,21,22,23,24];
ques_opt_incl.alcohol   = [];
ques_opt_incl.apathy    = [1,2,7,8,16,17,18];
ques_opt_incl.eat       = [1,11];
ques_opt_incl.anxiety   = [1,3,4,5,7,8,10,11,12,13,16,17,19,20];
ques_opt_incl.icar      = 1:16;

% count number of subjects with weird number of responses
subjctr = 0;
for i = 1:nsubj
    c = cell(size(test,1),1);
    c(:) = {i};
    total = sum(cell2mat(cellfun(@eq,test(:,1),c,'UniformOutput',false)));
    if total ~= nques
        if total > 0
            fprintf('subject %d has %d/%d questionnaires\n',i,total,nques);
        end
    elseif total >= nques
        subjctr = subjctr+1;
    end
end

ques_table = test;

% check that subjects with normal number of responses are complete
isubj = 0;
nrow  = size(ques_table,1);
isubj_incompl = [];
label_incompl = {};

ctr_incompl = 0;
for irow = 1:nrow
    if ques_table{irow,1} ~= isubj
        isubj = ques_table{irow,1};
    else
        continue
    end
    c    = cell(nrow,1);
    c(:) = {isubj};
    ind_subj = cell2mat(cellfun(@eq,ques_table(:,1),c,'UniformOutput',false));
    idx_row = find(ind_subj==1);
    
    labelctr = zeros(size(labelcnt));
    % adapted for RLINF code
    for ilabel = 1:numel(labelstr)
        labelctr(ilabel) = sum(contains(ques_table(idx_row,2),labelstr{ilabel}));
    end

    missedlabel = ~(labelctr == labelcnt);
    if any(missedlabel)
        fprintf('Subject %d does not have a full questionnaire dataset!\n',isubj);
        isubj_incompl = [isubj_incompl isubj];
        ctr_incompl = ctr_incompl + 1;
        missedlabelstr = {};
        for ilabel = 1:numel(labelctr)
            if missedlabel(ilabel)
                 missedlabelstr{end+1} =  labelstr{ilabel};
            end
        end
        label_incompl{ctr_incompl} = missedlabelstr;
    end
end
fprintf('\nAll questionnaires accounted for.\n');
ques_struct = cell(nsubj,1);

%% Process questionnaire: ALCOHOL
% find indices in ques_table with the chosen questionnaire
label = 'alcohol';
if ~is_short_ques
    nresp = 10;
    ind_test = ismember(ques_table(:,2),{label});
    idx_test = find(ind_test == 1)';
    subjlist_alc = [];
    
    % find excluded subjects (for the specific questionnaire)
    excl     = [];
    for iinc = 1:numel(label_incompl)
        if any(ismember(label_incompl{iinc},{label}))
            excl = [excl isubj_incompl(iinc)];
        end
    end
    
    for idx = idx_test
        isubj = ques_table{idx,1};
        if ismember(isubj,excl)
            continue
        end
        subjlist_alc = [subjlist_alc; isubj];
        if isempty(ques_struct{isubj})
            ques_struct{isubj} = struct;
        end
        ques_struct{isubj}.alcohol = struct;
        ques_struct{isubj}.alcohol.raw  = [];
        resp_str = ques_table(idx,3:nresp+2);
        [~,resp_str] = strtok(resp_str,':"');
        for iresp = 1:nresp
            if str2double(resp_str{iresp}(3)) == 0
                fprintf('Subject %d has missing questionnaire (''%s'')\n',isubj,label);
                isubj_incompl = [isubj_incompl isubj];
                continue_flag = true;
            break
        end
            ques_struct{isubj}.alcohol.raw = [ques_struct{isubj}.alcohol.raw str2double(resp_str{iresp}(3))];
        end
        
        % ALCOHOL SCORING
        %   Questions 1-8 are scored from 0 to 4
        %   Questions 9 and 10 are scored 0, 2, or 4
        alcohol_a = ques_struct{isubj}.alcohol.raw;
        alcohol_a = alcohol_a-1; % rescale from [1,5] to [0,4]
        
        for y = 9:10
            if alcohol_a(y) == 1
                alcohol_a(y) = 2;
            elseif alcohol_a(y) == 2
                alcohol_a(y) = 4;
            end
        end
        alcohol_s = sum(alcohol_a);
        
        ques_struct{isubj}.alcohol.raw   = alcohol_a;
        ques_struct{isubj}.alcohol.score = alcohol_s;
    end
else
    fprintf('Questionnaire (%s) excluded from analysis..\n',label);
end
fprintf('Processed questionnaire (%s)...\n',label);

%% Process questionnaire: ANXIETY
% find indices in ques_table with the chosen questionnaire
label = 'anxiety';
nresp_full = 20;
if ~is_short_ques
    nresp = nresp_full;
else
    nresp = numel(ques_opt_incl.(label));
end

ind_test = ismember(ques_table(:,2),{label});
idx_test = find(ind_test == 1)';
subjlist_anx = [];

% find excluded subjects (for the specific questionnaire)
excl     = [];
for iinc = 1:numel(label_incompl)
    if any(ismember(label_incompl{iinc},{label}))
        excl = [excl isubj_incompl(iinc)];
    end
end

for idx = idx_test
    isubj = ques_table{idx,1};
    if ismember(isubj,excl)
        continue
    end
    subjlist_anx = [subjlist_anx; isubj];
    if isempty(ques_struct{isubj})
        ques_struct{isubj} = struct;
    end
    ques_struct{isubj}.anxiety      = struct;
    ques_struct{isubj}.anxiety.raw  = nan(1,nresp_full);
    
    resp_str = ques_table(idx,3:nresp+2);
    [~,resp_str] = strtok(resp_str,':"');
    for iresp = 1:nresp
        if str2double(resp_str{iresp}(3)) == 0
            fprintf('Subject %d has missing questionnaire (''%s'')\n',isubj,label);
            isubj_incompl = [isubj_incompl isubj];
            continue_flag = true;
            break
        end
        if ~is_short_ques
            ques_struct{isubj}.anxiety.raw(iresp) = str2double(resp_str{iresp}(3));
        else
            ques_struct{isubj}.anxiety.raw(ques_opt_incl.(label)(iresp)) = str2double(resp_str{iresp}(3));
        end
    end
    
    % ANXIETY SCORING
    anxiety_a           = ques_struct{isubj}.anxiety.raw;
    anxiety_revitems    = [1 3 6 7 10 13 14 16 19];
    anxiety_qns         = zeros(1,20);
    for anxiety_rev = 1:length(anxiety_revitems)
        anxiety_qns(anxiety_revitems(anxiety_rev))  = 5 - anxiety_a(anxiety_revitems(anxiety_rev));
        anxiety_a(anxiety_revitems(anxiety_rev))    = 0;
    end
    anxiety_s = nansum(anxiety_a + anxiety_qns);
    
    ques_struct{isubj}.anxiety.raw   = anxiety_a + anxiety_qns;
    ques_struct{isubj}.anxiety.score = anxiety_s;
    if is_short_ques
        ques_struct{isubj}.anxiety.description = 'Score calculated from shortened/optimized questionnaire.';
    end
end
fprintf('Processed questionnaire (%s)...\n',label);
%% Process questionnaire: APATHY
% find indices in ques_table with the chosen questionnaire
label = 'apathy';
nresp_full = 18;
if ~is_short_ques
    nresp = nresp_full;
else
    nresp = numel(ques_opt_incl.(label));
end

ind_test = ismember(ques_table(:,2),{label});
idx_test = find(ind_test == 1)';
subjlist_apa = [];

% find excluded subjects (for the specific questionnaire)
excl     = [];
for iinc = 1:numel(label_incompl)
    if any(ismember(label_incompl{iinc},{label}))
        excl = [excl isubj_incompl(iinc)];
    end
end

for idx = idx_test
    isubj = ques_table{idx,1};
    if ismember(isubj,excl)
        continue
    end
    subjlist_apa = [subjlist_apa; isubj];
    if isempty(ques_struct{isubj})
        ques_struct{isubj}              = struct;
    end
    ques_struct{isubj}.apathy      = struct;
    ques_struct{isubj}.apathy.raw  = nan(1,nresp_full);
    
    resp_str = ques_table(idx,3:nresp+2);
    [~,resp_str] = strtok(resp_str,':"');
    
    for iresp = 1:nresp
        if str2double(resp_str{iresp}(3)) == 0
            fprintf('Subject %d has missing questionnaire (''%s'')\n',isubj,label);
            isubj_incompl = [isubj_incompl isubj];
            continue_flag = true;
            break
        end
        if ~is_short_ques
            ques_struct{isubj}.apathy.raw(iresp) = str2double(resp_str{iresp}(3));
        else
            ques_struct{isubj}.apathy.raw(ques_opt_incl.(label)(iresp)) = str2double(resp_str{iresp}(3));
        end
    end
    
    % APATHY SCORING
    %   apathy is coded so higher number means higher apathy, so reverse score the
    %   entire thing except for qn 6, 10 and 11
    apathy_a        = ques_struct{isubj}.apathy.raw;
    apathy_revitems = [1 2 3 4 5 7 8 9 12 13 14 15 16 17 18];
    apathy_qns      = zeros(1,18);
    for apathy_rev = 1:length(apathy_revitems)
        apathy_qns(apathy_revitems(apathy_rev)) = 5 - apathy_a(apathy_revitems(apathy_rev));
        apathy_a(apathy_revitems(apathy_rev))   = 0;
    end
    apathy_s = nansum(apathy_a+apathy_qns);
    
    ques_struct{isubj}.apathy.raw   = apathy_a+apathy_qns;
    ques_struct{isubj}.apathy.score = apathy_s;
    if is_short_ques
        ques_struct{isubj}.apathy.description = 'Score calculated from shortened/optimized questionnaire.';
    end
end
fprintf('Processed questionnaire (%s)...\n',label);

%% Process questionnaire: BIS
% find indices in ques_table with the chosen questionnaire
label = 'bis';
nresp_full = 30;
if ~is_short_ques
    nresp = nresp_full;
else
    nresp = numel(ques_opt_incl.(label));
end

ind_test = ismember(ques_table(:,2),{label});
idx_test = find(ind_test == 1)';
subjlist_bis = [];

% find excluded subjects (for the specific questionnaire)
excl     = [];
for iinc = 1:numel(label_incompl)
    if any(ismember(label_incompl{iinc},{label}))
        excl = [excl isubj_incompl(iinc)];
    end
end

for idx = idx_test
    isubj = ques_table{idx,1};
    if ismember(isubj,excl)
        continue
    end
    subjlist_bis = [subjlist_bis; isubj];
    if isempty(ques_struct{isubj})
        ques_struct{isubj}          = struct;
    end
    ques_struct{isubj}.bis      = struct;
    ques_struct{isubj}.bis.raw  = nan(1,nresp_full);
        
    resp_str = ques_table(idx,3:nresp+2);
    [~,resp_str] = strtok(resp_str,':"');
    for iresp = 1:nresp
        if str2double(resp_str{iresp}(3)) == 0
            fprintf('Subject %d has missing questionnaire (''%s'')\n',isubj,label);
            isubj_incompl = [isubj_incompl isubj];
            continue_flag = true;
            break
        end
        if ~is_short_ques
            ques_struct{isubj}.bis.raw(iresp) = str2double(resp_str{iresp}(3));
        else
            ques_struct{isubj}.bis.raw(ques_opt_incl.(label)(iresp)) = str2double(resp_str{iresp}(3));
        end
    end
    
    % BIS SCORING
    bis_a       = ques_struct{isubj}.bis.raw;
    bisrevitems = [1 7 8 9 10 12 13 15 20 29 30];
    for bis_rev = 1: length(bisrevitems)
        bis_qns(bisrevitems(bis_rev)) = 5 - bis_a(bisrevitems(bis_rev));
        bis_a(bisrevitems(bis_rev)) = 0;
    end
    bis_total = bis_qns+bis_a';
    
    % first order
    bis_atten       = [11 28 5 9 20];
    bis_atten_s     = nansum(bis_total(bis_atten));
    bis_motrimpul	= [17 19 22 3 2 25 4];
    bis_motrimpul_s = nansum(bis_total(bis_motrimpul));
    bis_selftrl     = [12 1 8 7 13 14];
    bis_selftrl_s   = nansum(bis_total(bis_selftrl));
    bis_cogcompl    = [15 29 10 27 18];
    bis_cogcompl_s  = nansum(bis_total(bis_cogcompl));
    bis_persev      = [21 16 30 23];
    bis_persev_s    = nansum(bis_total(bis_persev));
    bis_coginsta    = [26 6 24];
    bis_coginsta_s  = nansum(bis_total(bis_coginsta));
    % second order
    bis_attenimpul      = bis_atten_s + bis_coginsta_s;
    bis_motorimpul      = bis_motrimpul_s + bis_persev_s;
    bis_nonplanimpul    = bis_selftrl_s + bis_cogcompl_s;
    
    bis_s = bis_attenimpul + bis_motorimpul + bis_nonplanimpul;
    
    ques_struct{isubj}.bis.score.total                  = bis_s;
    ques_struct{isubj}.bis.score.firstord.atten         = bis_atten_s;
    ques_struct{isubj}.bis.score.firstord.mtrimpul      = bis_motrimpul_s;
    ques_struct{isubj}.bis.score.firstord.selfctrl      = bis_selftrl_s;
    ques_struct{isubj}.bis.score.firstord.cogcompl      = bis_cogcompl_s ;
    ques_struct{isubj}.bis.score.firstord.persev        = bis_persev_s;
    ques_struct{isubj}.bis.score.firstord.coginsta      = bis_coginsta_s;
    ques_struct{isubj}.bis.score.secondord.attenimpul   = bis_attenimpul;
    ques_struct{isubj}.bis.score.secondord.motorimpul   = bis_motorimpul;
    ques_struct{isubj}.bis.score.secondord.nonplanimpul = bis_nonplanimpul;
    if is_short_ques
        ques_struct{isubj}.bis.description = 'Score of all order calculated from shortened/optimized questionnaire.';
    end
end
fprintf('Processed questionnaire (%s)...\n',label);

%% Process questionnaire: DEPRESS
% find indices in ques_table with the chosen questionnaire
label = 'depress';
nresp_full = 20;
if ~is_short_ques
    nresp = nresp_full;
else
    nresp = numel(ques_opt_incl.(label));
end

ind_test = ismember(ques_table(:,2),{label});
idx_test = find(ind_test == 1)';
subjlist_dep = [];

% find excluded subjects (for the specific questionnaire)
excl     = [];
for iinc = 1:numel(label_incompl)
    if any(ismember(label_incompl{iinc},{label}))
        excl = [excl isubj_incompl(iinc)];
    end
end

for idx = idx_test
    isubj = ques_table{idx,1};
    if ismember(isubj,excl)
        continue
    end
    subjlist_dep = [subjlist_dep; isubj];
    if isempty(ques_struct{isubj})
        ques_struct{isubj} = struct;
    end
    ques_struct{isubj}.depress     = struct;
    ques_struct{isubj}.depress.raw = nan(1,nresp_full);
    
    resp_str = ques_table(idx,3:nresp+2);
    [~,resp_str] = strtok(resp_str,':"');
    for iresp = 1:nresp
        if str2double(resp_str{iresp}(3)) == 0
            fprintf('Subject %d has missing questionnaire (''%s'')\n',isubj,label);
            isubj_incompl = [isubj_incompl isubj];
            continue_flag = true;
            break
        end
        if ~is_short_ques
            ques_struct{isubj}.depress.raw(iresp) = str2double(resp_str{iresp}(3));
        else
            ques_struct{isubj}.depress.raw(ques_opt_incl.(label)(iresp)) = str2double(resp_str{iresp}(3));
        end
    end
    
    % DEPRESS SCORING
    depress_a       = ques_struct{isubj}.depress.raw;
    zung_revitems   = [2 5 6 11 12 14 16 17 18 20];
    for zung_rev = 1: length(zung_revitems)
        depress_qns(zung_revitems(zung_rev)) = 5 - depress_a(zung_revitems(zung_rev));
        depress_a(zung_revitems(zung_rev)) = 0;
    end
    depress_s = nansum(depress_a+depress_qns);
    
    ques_struct{isubj}.depress.raw   = depress_a+depress_qns;
    ques_struct{isubj}.depress.score = depress_s;
    if is_short_ques
        ques_struct{isubj}.depress.description = 'Score calculated from shortened/optimized questionnaire.';
    end
end
fprintf('Processed questionnaire (%s)...\n',label);

%% Process questionnaire: EAT
% find indices in ques_table with the chosen questionnaire
label = 'eat';
nresp_full = 26;
if ~is_short_ques
    nresp = nresp_full;
else
    nresp = numel(ques_opt_incl.(label));
end

ind_test = ismember(ques_table(:,2),{label});
idx_test = find(ind_test == 1)';
subjlist_eat = [];

% find excluded subjects (for the specific questionnaire)
excl     = [];
for iinc = 1:numel(label_incompl)
    if any(ismember(label_incompl{iinc},{label}))
        excl = [excl isubj_incompl(iinc)];
    end
end

for idx = idx_test
    isubj = ques_table{idx,1};
    if ismember(isubj,excl)
        continue
    end
    subjlist_eat = [subjlist_eat; isubj];
    if isempty(ques_struct{isubj})
        ques_struct{isubj} = struct;
    end
    ques_struct{isubj}.eat      = struct;
    ques_struct{isubj}.eat.raw  = nan(1,nresp_full);
    
    resp_str = ques_table(idx,3:nresp+2);
    [~,resp_str] = strtok(resp_str,':"');
    for iresp = 1:nresp
        if str2double(resp_str{iresp}(3)) == 0
            fprintf('Subject %d has missing questionnaire (''%s'')\n',isubj,label);
            isubj_incompl = [isubj_incompl isubj];
            continue_flag = true;
            break
        end
        if ~is_short_ques
            ques_struct{isubj}.eat.raw(iresp) = str2double(resp_str{iresp}(3));
        else
            ques_struct{isubj}.eat.raw(ques_opt_incl.(label)(iresp)) = str2double(resp_str{iresp}(3));
        end
    end
    
    % EAT SCORING
    eat_a = ques_struct{isubj}.eat.raw;
    
    eat_qns         = zeros(1,26);
    eat_qns(25)     = eat_a(25)-3;
    eat_a           = 4 - eat_a;
    eat_a(25)       = 0;
    eat_a           = eat_a+eat_qns;
    eat_a(eat_a<0)  = 0;

    eat_s = sum(eat_a);
    
    eat_diet        = [1 6 7 10 11 12 14 16 17 22 23 24 25];
    eat_diet_s      = nansum(eat_a(eat_diet));
    eat_bulimia     = [3 4 9 18 21 26];
    eat_bulimia_s   = nansum(eat_a(eat_bulimia));
    eat_oral        = [2 5 8 13 15 19 20];
    eat_oral_s      = nansum(eat_a(eat_oral));
    
    ques_struct{isubj}.eat.raw              = eat_a;
    ques_struct{isubj}.eat.score.total      = eat_s;
    ques_struct{isubj}.eat.score.oral       = eat_oral_s;
    ques_struct{isubj}.eat.score.bulimia    = eat_bulimia_s;
    ques_struct{isubj}.eat.score.diet       = eat_diet_s;
    if is_short_ques
        ques_struct{isubj}.eat.description = 'Score of all order calculated from shortened/optimized questionnaire.';
    end
end
fprintf('Processed questionnaire (%s)...\n',label);

%% Process questionnaire: IQ (will bug if trying to run twice w/o clearing)
% find indices in ques_table with the chosen questionnaire
label = 'icar';
ind_test = contains(ques_table(:,2),{label});
subjlist_iq = [];

% find excluded subjects (for the specific questionnaire)
excl = [];
% check for missing data
isone = true;
istwo = false;
isthr = false;
for i = 1:numel(ind_test)
    if istwo
        istwo = false;
        isthr = true;
        continue
    elseif isthr
        isthr = false;
        isone = true;
        continue
    end
    if ind_test(i) == 1
        arraythree = ind_test(i:i+2);
        if sum(arraythree) ~= 3
            excl = [excl; ques_table{i,1}];
        end
        isone = false;
        istwo = true;
    end
end

% extract iq scores
idx_test    = find(ind_test == 1)';
for idx = idx_test
    isubj = ques_table{idx,1};
    if ismember(isubj,excl)
        continue
    end
    subjlist_iq = union(subjlist_iq,isubj);
    if isempty(ques_struct{isubj})
        ques_struct{isubj} = struct;
    end
    if ~isfield(ques_struct{isubj},'iq')
        ques_struct{isubj}.iq       = struct;
        ques_struct{isubj}.iq.raw   = [];
    end
    
    if ~isempty(ques_table{idx,10}) % 7 responses
        resp_str = ques_table(idx,3:10);
        [~,resp_str] = strtok(resp_str,':"');
        for jq = 1:8
            ques_struct{isubj}.iq.raw = [ques_struct{isubj}.iq.raw str2double(resp_str{jq}(3))];
        end 
    else % 4 responses
        resp_str = ques_table(idx,3:6);
        [~,resp_str] = strtok(resp_str,':"');
        for jq = 1:4
            ques_struct{isubj}.iq.raw = [ques_struct{isubj}.iq.raw str2double(resp_str{jq}(3))];
        end 
    end
    if is_short_ques
        ques_struct{isubj}.iq.description = 'ICAR from the end of shortened questionnaire.';
    end
end

% IQ SCORING
iq_ans = [4 4 4 6 6 3 4 4 5 2 2 4 3 2 6 7];
iqs = nan(nsubj,1);
for isubj = subjlist_iq
    ques_struct{isubj}.iq.score = sum(bsxfun(@eq,ques_struct{isubj}.iq.raw,iq_ans));
    iqs(isubj) = ques_struct{isubj}.iq.score;
end
fprintf('Processed questionnaire (%s)...\n',label);

% save('icar_score.mat','iqs');

%% Process questionnaire: OCIR
% find indices in ques_table with the chosen questionnaire
label = 'ocir';
nresp_full = 19;
if ~is_short_ques
    nresp = nresp_full;
else
    nresp = numel(ques_opt_incl.(label));
end

ind_test = ismember(ques_table(:,2),{label});
idx_test = find(ind_test == 1)';
subjlist_ocir = [];
excl_catch_ocir = [];

% find excluded subjects (for the specific questionnaire)
excl     = [];
for iinc = 1:numel(label_incompl)
    if any(ismember(label_incompl{iinc},{label}))
        excl = [excl isubj_incompl(iinc)];
    end
end

for idx = idx_test
    isubj = ques_table{idx,1};
    if ismember(isubj,excl)
        continue
    end
    continue_flag = false; % to skip analysis of subjects with faulty/inattentive data
    subjlist_ocir = [subjlist_ocir; isubj];
    if isempty(ques_struct{isubj})
        ques_struct{isubj} = struct;
    end
    ques_struct{isubj}.ocir      = struct;
    ques_struct{isubj}.ocir.raw  = nan(1,nresp_full);
    
    resp_str = ques_table(idx,3:nresp+2);
    [~,resp_str] = strtok(resp_str,':"');
    for iresp = 1:nresp
        % check for submitted questionnaires with no data
        if str2double(resp_str{iresp}(3)) == 0
            fprintf('Subject %d has missing questionnaire (''%s'')\n',isubj,label);
            isubj_incompl = [isubj_incompl isubj];
            continue_flag = true;
            break
        end
        if ~is_short_ques
            ques_struct{isubj}.ocir.raw(iresp) = str2double(resp_str{iresp}(3));
        else
            ques_struct{isubj}.ocir.raw(ques_opt_incl.(label)(iresp)) = str2double(resp_str{iresp}(3));
        end
    end
    if continue_flag
        continue
    end
    
    % OCIR SCORING
    ocir_a  = ques_struct{isubj}.ocir.raw;
    catchqn = ocir_a([12]); % number 12 is the catch question, maybe save this data too later
    ocir_a  = ocir_a([1:11 13:19]);
    if catchqn ~= 2
        fprintf('Subject %d missed catch item\n',isubj);
        excl_catch_ocir = [excl_catch_ocir isubj];
        isubj_incompl = [isubj_incompl isubj];
    end
    ocir_qns    = (ocir_a-1);
    ocir_s      = nansum(ocir_qns);
    
    ques_struct{isubj}.ocir.raw   = ocir_qns;
    ques_struct{isubj}.ocir.score = ocir_s;
    if is_short_ques
        ques_struct{isubj}.ocir.description = 'Score calculated from shortened/optimized questionnaire.';
    end
end
fprintf('Processed questionnaire (%s)...\n',label);

%% Process questionnaire: SCHIZO
% find indices in ques_table with the chosen questionnaire
label = 'schizo';
if ~is_short_ques   
    ind_test = contains(ques_table(:,2),{label});
    subjlist_sch = [];
    
    % find excluded subjects (for the specific questionnaire)
    excl = [];
    % check for missing data
    isone = true;
    istwo = false;
    for i = 1:numel(ind_test)
        if istwo
            istwo = false;
            isthr = true;
            continue
        end
        if ind_test(i) == 1
            arraytwo = ind_test(i:i+1);
            if sum(arraytwo) ~= 2
                excl = [excl; ques_table{i,1}];
            end
            isone = false;
            istwo = true;
        end
    end
    
    idx_test    = find(ind_test == 1)';
    for idx = idx_test
        isubj = ques_table{idx,1};
        if ismember(isubj,excl)
            continue
        end
        subjlist_sch = union(subjlist_sch,isubj);
        if isempty(ques_struct{isubj})
            ques_struct{isubj} = struct;
        end
        if ~isfield(ques_struct{isubj},'schizo')
            ques_struct{isubj}.schizo       = struct;
            ques_struct{isubj}.schizo.raw   = [];
        end
        
        if ~isempty(ques_table{idx,25}) % 23 responses
            resp_str = ques_table(idx,3:25);
            [~,resp_str] = strtok(resp_str,':"');
            for jq = 1:23
                ques_struct{isubj}.schizo.raw = [ques_struct{isubj}.schizo.raw str2double(resp_str{jq}(3))];
            end 
        else % 20 responses
            resp_str = ques_table(idx,3:22);
            [~,resp_str] = strtok(resp_str,':"');
            for jq = 1:20
                ques_struct{isubj}.schizo.raw = [ques_struct{isubj}.schizo.raw str2double(resp_str{jq}(3))];
            end 
        end
    end
    
    % SCHIZO SCORING
    schizo_qns1     = zeros(1,43);
    schizo_revitems = [26 27 28 30 31 34 37 39];
    for isubj = subjlist_sch
        schizo_a    = ques_struct{isubj}.schizo.raw;
        schizo_qns  = schizo_a - 1;
        for schizo_rev = 1: length(schizo_revitems) 
            schizo_qns1(schizo_revitems(schizo_rev)) = 1 - schizo_qns(schizo_revitems(schizo_rev));
            schizo_qns(schizo_revitems(schizo_rev)) = 0;
        end
        schizo_total        = schizo_qns1 + schizo_qns;
        schizo_unuslexp     = sum(schizo_total(1:12));
        schizo_cogdisorg    = sum(schizo_total(13:23));
        schizo_introanhed   = sum(schizo_total(24:33));
        schizo_impulnoncon  = sum(schizo_total(34:43));
    
        schizo_s = sum(schizo_total);
    
        ques_struct{isubj}.schizo.raw               = schizo_total;
        ques_struct{isubj}.schizo.score.total       = schizo_s;
        ques_struct{isubj}.schizo.score.unuslexp    = schizo_unuslexp;
        ques_struct{isubj}.schizo.score.cogdisorg   = schizo_cogdisorg;
        ques_struct{isubj}.schizo.score.introanhed  = schizo_introanhed;
        ques_struct{isubj}.schizo.score.impulnoncon = schizo_impulnoncon;
    end
else
    fprintf('Questionnaire (%s) excluded from analysis..\n',label);
end
fprintf('Processed questionnaire (%s)...\n',label);

%% Process questionnaire: SOCIAL
% find indices in ques_table with the chosen questionnaire
label = 'social';
nresp_full = 48;
if ~is_short_ques
    nresp = nresp_full;
else
    nresp = numel(ques_opt_incl.(label));
    odd_idx = 1:2:nresp_full;
    eve_idx = 2:2:nresp_full;
end

ind_test = contains(ques_table(:,2),{label});
subjlist_soc = [];

% find excluded subjects (for the specific questionnaire)
excl = [];
% check for missing data
isone = true;
istwo = false;
for i = 1:numel(ind_test)
    if istwo
        istwo = false;
        isthr = true;
        continue
    end
    if ind_test(i) == 1
        arraytwo = ind_test(i:i+1);
        if sum(arraytwo) ~= 2
            excl = [excl; ques_table{i,1}];
        end
        isone = false;
        istwo = true;
    end
end

% filter and register item scores
idx_test    = find(ind_test == 1)';
isubj_prev = 0; % used for shortened questionnaire filtering
for idx = idx_test
    isubj = ques_table{idx,1};
    if ismember(isubj,excl)
        continue
    end
    subjlist_soc = union(subjlist_soc,isubj);
    if isempty(ques_struct{isubj})
        ques_struct{isubj} = struct;
    end

    if ~isfield(ques_struct{isubj},'social')
        ques_struct{isubj}.social = struct;
    end
    
    if ~is_short_ques
        ques_struct{isubj}.social.raw = [];
        resp_str = ques_table(idx,3:26); % 24 responses
        [~,resp_str] = strtok(resp_str,':"');
        for jq = 1:24
            ques_struct{isubj}.social.raw = [ques_struct{isubj}.social.raw str2double(resp_str{jq}(3))];
        end 
    else
        if isubj == isubj_prev % 2nd page of social
            n_end = 22; % 10 items, 20 responses
            ptr_add = 9;
        else
            ques_struct{isubj}.social.raw = nan(1,nresp_full);
            n_end = 20; % 9 items, 18 responses
            ptr_add = 0;
        end
        resp_str = ques_table(idx,3:n_end);
        [~,resp_str] = strtok(resp_str,':"');

        for iresp = 1:n_end-2
            ptr = round(iresp/2)+ptr_add;
            if mod(iresp,2) == 1 % odd
                ques_struct{isubj}.social.raw(odd_idx(ques_opt_incl.(label)(ptr))) = str2double(resp_str{iresp}(3));
            else % even
                ques_struct{isubj}.social.raw(eve_idx(ques_opt_incl.(label)(ptr))) = str2double(resp_str{iresp}(3));
            end
        end
    end
    isubj_prev = isubj;
end

% SOCIAL SCORING
for isubj = subjlist_soc
    leb_a   = ques_struct{isubj}.social.raw;
    leb_qns = leb_a-1;
    leb_s   = sum(leb_qns);

    y = 1;
    x = 2;
    for z = 1:24
        leb_qns_avg(z) = (leb_qns(y)+leb_qns(x))/2;
        y = y + 2;
        x = x + 2;
    end

    ques_struct{isubj}.social.raw	= leb_qns;
    ques_struct{isubj}.social.avg   = leb_qns_avg;
    ques_struct{isubj}.social.score = leb_s;
end
fprintf('Processed questionnaire (%s)...\n',label);

%% Organize

% exclude subjects w/o a full questionnaire and who failed the catch
% question in OCIR
idx_excl_ques = union(isubj_incompl,excl_catch_ocir);
fprintf('\nQuestionnaire-excluded subjects:\n');
fprintf(' %d,',idx_excl_ques);
fprintf('\n');

save(sprintf('./processed/%s/idx_excl_ques.mat',samplename),'idx_excl_ques');
save(sprintf('./processed/%s/ques_struct.mat',samplename),'ques_struct'); % save the raw questionnaire scores

fprintf('Questionnaire data saved.\n');

