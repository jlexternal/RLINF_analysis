% test_simulate_kf
clear all

samplename = 'pilot07'; % Pilot - 'pilot'

load(sprintf('./processed/%s/preprocessed_data_%s.mat',samplename,samplename)); % load the raw data structure for data sample

isubj = 1;

cfg = struct;
fbabs = idx_fbabs(isubj,:);
fb    = idx_fb(isubj,:);
bmrsp = idx_blmn(1,:);
bmst  = idx_bmstate(isubj,:);

% Need to create 2 vectors for the values of options 1 and 2
% blue/moon is option 1
fb_c = fbabs';
fb_i = 100-fb_c;
fb_ci = [fb_c fb_i];

nt = size(fb_ci,1);
fb1 = fb_ci(sub2ind(size(fb_ci),(1:nt)',-bmst'+2));
fb_opt = [fb1 100-fb1]';

% filter by condition
for icond = 0:1
    idx = idx_cond(isubj,:) == icond;

    trl = idx_trial(idx);
    blk = idx_blk(idx);
    r1 = fb_opt(1,idx);
    r2 = fb_opt(2,idx);
end

simulate_kf
