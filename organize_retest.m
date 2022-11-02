% organize_retest
clear all

samplename = 'sample2';
versionadd = '';
dirname = './import';
filename = sprintf('pid_%s%s.csv',samplename,versionadd);

dir = ls(fullfile(sprintf('%s/%s/%s', dirname, samplename, filename)));

fulltable = readtable(dir(1:end-1));
dat = load(sprintf('./processed/%s/preprocessed_data_%s.mat',samplename,samplename));
load(sprintf('./processed/%s/idx_excl_ques.mat',samplename),'idx_excl_ques');
load(sprintf('./processed/%s/pcor_%s.mat',samplename,samplename),'pcorr');

idx_ques = ones(size(dat.idx_fb(:,1),1),1);
idx_ques(idx_excl_ques') = 0;
idx_task = ~isnan(dat.idx_fb(:,1));

p_threshold = 158/288; % 158 out of 288 for binomial test against random

idx_aboveChance = any(pcorr >= p_threshold,2);

idx_retest = idx_task & idx_ques & idx_aboveChance;
