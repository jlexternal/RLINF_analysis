% gen_task_testcases

% Generates a task with stereotyped settings to produce predictable results
% from model outputs.

% Test case blocks:
% 1/ Stationary values at extremes
% 2/ Stationary values at midpoint (e.g. 50)
% 3/ Stationary values at moderate values
% 4/ Noisy values at moderate values
% 5/ Volatile (slow change-point) values with no ambiguity
% 6/ Volatile (rapid change-point) values with no ambiguity
% 7/ Volatile (slow drifting) values with no ambiguity
% 8/ Volatile (rapid drifting) values with no ambiguity

% -------------- Input --------------
% Test parameters:
ncases = 8; % number of test cases to test
ntrl_per_state  = 12; % number of trials for a state 
savename = 'task_testcases_rlinf';
% -----------------------------------

testcases = cell(ncases,1);
% Task generator loop
for icase = 1:ncases
    isvolatile = false;
    isdrifting = false;
    israpid    = false;
    isnoisy    = false;

    % set configuration
    switch icase
        case 1 % extreme value stationary blocks
            descript = 'Stationary values at extremes';
            casevals = [0 1];
        case 2 % midpoint value stationary blocks
            descript = 'Stationary values at midpoint (e.g. 50)';
            casevals = .5;
        case 3 % moderate value stationary blocks
            descript = 'Stationary values at moderate values';
            casevals = [.25 .4 .6 .75];
        case 4 % moderate value noisy blocks
            descript = 'Noisy values at moderate values';
            casevals = [.25 .4 .6 .75];
            isnoisy  = true;
        case 5 % moderate value stationary slow volatile blocks 
            descript = 'Volatile (slow change-point) values with no ambiguity';
            casevals = [.25 .75 .3 .6];
            isvolatile = true;
        case 6 % moderate value stationary rapid volatile blocks
            descript = 'Volatile (rapid change-point) values with no ambiguity';
            casevals = [.25 .75 .3 .6];
            isvolatile = true;
            israpid = true;
        case 7 % slow drifting value volatile blocks
            descript = 'Volatile (slow drifting) values with no ambiguity';
            casevals = [.25 .75 .3 .6];
            isvolatile = true;
            isdrifting = true;
        case 8 % rapid drifting value volatile blocks
            descript = 'Volatile (rapid drifting) values with no ambiguity';
            casevals = [.25 .75 .3 .6];
            isvolatile = true;
            isdrifting = true;
            israpid = true;
        otherwise
            warning('Case %d is not accounted for in generator loop!',icase);
    end

    % Generate trials based on configuration above
    fprintf('Generating task values for description:\n> "%s"\n',descript);

    ncasevals = numel(casevals);
    if isvolatile
        if isdrifting
            idx_trl = 1:(ncasevals-1)*ntrl_per_state;
            idx_fb = nan(1,(ncasevals-1)*ntrl_per_state);
        else
            idx_trl = 1:ncasevals*ntrl_per_state;
            idx_fb = nan(1,ncasevals*ntrl_per_state);
        end
    else
        idx_trl = nan(1,ncasevals*ntrl_per_state);
        idx_fb = nan(1,ncasevals*ntrl_per_state);
    end
    
    ptr = 1; % trial index pointer
    for ival = 1:ncasevals
        val = casevals(ival); 
        idx = ptr:ptr+ntrl_per_state-1;
        
        if ~isvolatile % stationary cases (N = 4)
            idx_trl(idx) = 1:ntrl_per_state;
            idx_fb(idx)  = ones(ntrl_per_state,1)*val;
        else % volatile cases (N = 4)
            if isdrifting
                if val == casevals(end)
                    % volatile cases will look at 2 values (i and i+1) simulataneously
                    break
                end
                idx_fb(idx) = linspace(val,casevals(ival+1),ntrl_per_state);
            else
                idx_fb(idx) = val;
            end
        end
        ptr = ptr + ntrl_per_state;
    end
    if israpid
        idx_fb(mod(1:numel(idx_fb),2)==0) = [];
        idx_trl = 1:numel(idx_fb);
    end
    if isnoisy
        idx_fb = normrnd(idx_fb,.05);
    end

    testcases{icase} = struct;
    testcases{icase}.description = descript;
    testcases{icase}.idx_fb = idx_fb;
    testcases{icase}.idx_trl = idx_trl;
end
fprintf('All tasks generated (n = %d).\n',ncases);
save(sprintf('%s.mat',savename),'testcases');
fprintf('Test cases saved as ''%s.mat''\n',savename);
