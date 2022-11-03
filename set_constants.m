% set_constants
clear all
samplename = 'sample2';
ncnd = 2;
nblk = 4;
ntrl = 73;
nepis = 24;
fnr = .30;

condrgb = [93 74 25; 25 42 68]/100;

save(sprintf('./constants/constants_rlinf_%s',samplename));

% set indexes for task and questionnaires
load(sprintf('./processed/%s/preprocessed_data_%s.mat',samplename,samplename),'idx_blmn');
idx_fullTask = ~isnan(idx_blmn(:,1));

load(sprintf('./processed/%s/ques_struct.mat',samplename));
load(sprintf('./processed/%s/idx_excl_ques.mat',samplename));
idx_fullQues = ~cellfun(@isempty,ques_struct);
idx_fullQues(idx_excl_ques) = 0;

idx_fullAll = idx_fullTask & idx_fullQues;
save(sprintf('./processed/%s/idx_TaskQuesAll.mat',samplename),'idx_fullTask','idx_fullQues','idx_fullAll');