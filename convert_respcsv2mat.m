% convert_respcsv2mat.m
%
% Usage: Converts csv files taken from the MySQL database (with predefined
%           selection clauses) into MATLAB structures.
% 
% Note: This script only works for files generated from Pilot v6 and on.
clear all
clc

% -------------- Input --------------
isretest = true;
samplename = 'sample2';
versionadd = '';
% -----------------------------------
dirname = './import';
filename = sprintf('resp_%s%s.csv',samplename,versionadd);
if isretest
    filename = sprintf('resp_retest_%s%s.csv',samplename,versionadd);
    load(sprintf('./processed/%s/idx_TaskQuesAll.mat',samplename),'idx_fullAll');
    load(sprintf('./processed/%s/idx_retested.mat',samplename),'idx_retested');
end

dir = ls(fullfile(sprintf('%s/%s/%s', dirname, samplename, filename)));
fulltable = readtable(dir(1:end-1));

% Check subjects for missing data points
ntrl_max = 584;

extra_trial_idx_fairy = [1 73 145 217];
extra_trial_idx_bandit = [72 144 216 288];

nsubj = max(fulltable.subj);
isubj_task_excl = [];
isubj_task_extra = [];
subj_struct = cell(1,nsubj);
clc
for isubj = 1:nsubj
    idx = fulltable.subj == isubj;
    
    % if number of trials are mismatched
    if sum(idx) ~= ntrl_max
        if sum(idx) == 0 % returned subjects
            isubj_task_excl = [isubj_task_excl isubj];
            continue
        end
        fprintf('Subject %d reported %d/%d trials!\n',isubj,sum(idx),ntrl_max);
        if sum(idx) < ntrl_max % find which trials are missing for partial trial sums
            isubj_task_excl = [isubj_task_excl isubj];
            % fairy has extra trial at the BEGINNING
            t_missing = setdiff(1:288,fulltable.itrl(idx & fulltable.icnd == 1));
            if ~isempty(t_missing)
                fprintf(' > %d fairy trials missing\n',numel(t_missing));
                if numel(t_missing) < 8;
                    fprintf('   (%s)',strjoin(string(t_missing),', '));
                end
            end
            for it = extra_trial_idx_fairy
                if sum(idx & fulltable.icnd == 1 & fulltable.itrl == it) ~= 2
                    fprintf('Subject %d has trial missing on block %d in the fairy condition\n', isubj, find(extra_trial_idx_fairy==it));
                    isubj_task_excl = [isubj_task_excl isubj];
                end
            end
            
            % bandit has extra trial at the END
            t_missing = setdiff(1:288,fulltable.itrl(idx & fulltable.icnd == 0));
            if ~isempty(t_missing)
                fprintf(' > %d bandit trials missing\n',numel(t_missing));
                if numel(t_missing) < 8;
                    fprintf('   (%s)',strjoin(string(t_missing),', '));
                end
            end
            for it = extra_trial_idx_bandit
                if sum(idx & fulltable.icnd == 0 & fulltable.itrl == it) ~= 2
                    fprintf('Subject %d has trial missing on block %d in the bandit condition\n', isubj, find(extra_trial_idx_bandit==it));
                    isubj_task_excl = [isubj_task_excl isubj];
                end
            end
        else % subjects can have more than max number of trials
            warning('Subject %d has extra trials. Investigate...',isubj);
            % check for missing fairy trials
            t_missing = setdiff(1:288,fulltable.itrl(idx & fulltable.icnd == 1));
            if ~isempty(t_missing)
                fprintf(' > %d fairy trials missing\n',numel(t_missing));
                isubj_task_excl = [isubj_task_excl isubj];
            else
                fprintf(' > all fairy trials accounted for. Requires manual filtering.\n');
                isubj_task_extra = [isubj_task_extra isubj];
            end
            % check for missing bandit trials
            t_missing = setdiff(1:288,fulltable.itrl(idx & fulltable.icnd == 0));
            if ~isempty(t_missing)
                fprintf(' > %d bandit trials missing\n',numel(t_missing));
                isubj_task_excl = [isubj_task_excl isubj];
            else
                fprintf(' > all bandit trials accounted for. Requires manual filtering.\n');
                isubj_task_extra = [isubj_task_extra isubj];
            end
        end
        
        fprintf('\n');
        subj_struct{isubj} = fulltable(idx,:);
        continue
    end
    subj_struct{isubj} = fulltable(idx,:);
end

fprintf('Found %d/%d subjects with null/partial task dataset.\n',numel(unique(isubj_task_excl)),nsubj)

% remove data from subjects with partial or incomplete datasets
for isubj = unique(isubj_task_excl)
    subj_struct{isubj} = [];
end

% Filter extra trials if necessary
condstr = {'bandit','fairy'};
for isubj = unique(isubj_task_extra)
    isFixableDuplicate = true;
    while isFixableDuplicate
        dat = subj_struct{isubj};
        nt_subj = size(dat,1);
        fprintf('Analyzing data of subject %d...\n',isubj);
        % check which condition has extra trials
        for icond = 0:1
            nt = sum(dat.icnd == icond);
            if nt > ntrl_max/2
                fprintf('Extra trials (%d/%d) found in the %s (cond %d) condition!\n',nt,ntrl_max/2,condstr{icond+1},icond);
                cond2fix = icond;
            end
        end
        
        % check if the extra trials are due to repeated timestamp logs
        if sum(dat.ts(2:end) == dat.ts(1:end-1)) == (nt_subj-ntrl_max)
            fprintf(['Extra trials are strictly due to repeated timestamped entries!\n' ...
                'Running automated deletion of repeated entries...\n']);
            it2delete = [];
            for it = 2:nt_subj
                if dat.icnd(it) ~= cond2fix
                    continue
                end
                % find duplicate timestamp trials (2nd occurence)
                if dat.ts(it) == dat.ts(it-1)
                    it2delete = [it2delete it];
                end
            end
            % delete these entries
            dat(it2delete,:) = [];
            subj_struct{isubj} = dat;
            if size(dat,1) == ntrl_max
                fprintf('Subject %d data has been fully corrected, replacing entry in subj_struct...\n',isubj);
                isFixableDuplicate = false;
            end
        else
            warning('Subject %d has extra entries that are not due to timestamp duplicates! Check manually.',isubj);
            isFixableDuplicate = false;
        end
    end
    fprintf('Corrections to subject %d data complete.\n\n',isubj);
end

isubj_task_excl = unique(isubj_task_excl);
fprintf('%d full data points collected.\n\n',nsubj-numel(isubj_task_excl));

% Save to file
savenameadd = '';
if isretest
    savenameadd = 'retest';
    idx_task_rt = ~cellfun(@isempty,subj_struct);
    save(sprintf('./processed/%s/idx_task_rt.mat',samplename),"idx_task_rt");
end
savename = sprintf('subj_struct_%s_%s',savenameadd,samplename);
savedir = sprintf('./processed/%s', samplename);

if not(isfolder(savedir))
    mkdir(savedir)
end
save(sprintf('%s/%s',savedir, savename), 'subj_struct');
fprintf('File saved as %s.mat in %s\n',savename, savedir);