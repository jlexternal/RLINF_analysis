% organize_retest
clear all

samplename = 'sample2';
versionadd = '';
dirname = './import';
filename = sprintf('pid_%s%s.csv',samplename,versionadd);

dir = ls(fullfile(sprintf('%s/%s/%s', dirname, samplename, filename)));

fulltable = readtable(dir(1:end-1));
dat = load(sprintf('./processed/%s/preprocessed_data_%s.mat',samplename,samplename));
load(sprintf('./processed/%s/idx_TaskQuesAll.mat',samplename),'idx_fullAll');
load(sprintf('./processed/%s/pcor_%s.mat',samplename,samplename),'pcorr');

p_threshold = 158/288; % 158 out of 288 for binomial test against random

idx_aboveChance = any(pcorr >= p_threshold,2);

idx_retest = idx_fullAll & idx_aboveChance;

% table for database and appication
fulltable.retest = idx_retest;

basepay = 7.50;
%%
% text for Prolific (batch 1)
fprintf('PIDs for retest (Sample 2 - run 1)\n\n');
nbatch = 0;
for isubj = 1:110
    if idx_retest(isubj)
        fprintf('%s,\n',fulltable.pid{isubj});
        nbatch = nbatch + 1;
    end
end
fprintf('%d subjects in batch 1 \n\n',nbatch);
avgbonus = ((7*2) + (36))/nbatch;
fprintf('Participants gained an average bonus of %.2f in run 1.\n',avgbonus);

fprintf('Expect to pay %.2f for batch 1 retest...\n',nbatch*(basepay + avgbonus));

%%
% text for Prolific (batch 2)
fprintf('\nPIDs for retest (Sample 2 - run 2)\n\n');
nbatch = 0;

ctr = 0;
for isubj = 111:size(fulltable,1)
    if idx_retest(isubj)
        fprintf('%s,\n',fulltable.pid{isubj});
        nbatch = nbatch + 1;

        ctr = ctr + 1;
        if mod(ctr,25) == 0
            fprintf('\n')
        end

    end
end
fprintf('%d subjects in batch 1 \n\n',nbatch);
avgbonus = ((7*2) + (36))/nbatch;
fprintf('Participants gained an average bonus of %.2f in run 2.\n',avgbonus);

fprintf('Expect to pay %.2f for batch 2 retest...\n',nbatch*(basepay + avgbonus));


%% Generate SQL retest update query
clc
for isubj = 1:247
    if idx_retest(isubj)
        fprintf('UPDATE main_table SET retest = 1 WHERE subj = %d;\n',isubj);
    end
end

%% Get retest completion booleans

dirname = 'import';
samplename = 'sample2';
filename = 'pid_retest_sample2.csv';

dir = ls(fullfile(sprintf('./%s/%s/%s', dirname, samplename, filename)));
fulltable = readtable(dir(1:end-1));

idx_retested = fulltable.retest;
save(sprintf('./processed/%s/idx_retested.mat',samplename));


%%
% DROP PROCEDURE IF EXISTS countresp;
% DELIMITER $$
% CREATE PROCEDURE countresp(c1 INT, c2 INT)
% BEGIN
%     counter: LOOP
%         SET c1 = c1 + 1;
%         SET c2 = c2;
%         IF c1 < c2 THEN
%             SELECT COUNT(*) FROM resp_table WHERE uid LIKE CONCAT(LPAD(CONVERT(c1,CHAR),3,'0'),'%');
%             ITERATE counter;
%         END IF;
%         LEAVE counter;
%     END LOOP counter;
%     SET @x = c1;
% END $$
% DELIMITER ;