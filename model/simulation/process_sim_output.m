% process_sim_output
% 
% Function: Processes the output of simulations to facilitate easier
% analysis with pre-existing analysis scripts.
%
% Jun Seok Lee - October 2022

clear all

% specify conditions of simulation 
samplename = 'pilot07';
modeltype = 'noisyINF';
nsim = 1000;

filename = sprintf('out_sim_%s_nsim%d_%s',modeltype,nsim,samplename);

% load constants
load(sprintf('../../constants/constants_rlinf_%s.mat',samplename)); % load nblk, ncond, ntrl
% load simulation output
load(sprintf('./out/%s.mat',filename));

if ~strcmpi(sim_out.samplename,samplename)
    error('Sample ID of data does not match that of simulations!');
end
fprintf('\nProcessing file ''%s''...\n',filename);
fprintf('Processing simulations (n=%d) generated on %s (%s)...\n',nsim,datestr(sim_out.date),samplename);

sims = sim_out.out;
nsubj = size(sims,1);

isfirstpass = true;

for isubj = 1:nsubj
    if isempty(sims{isubj,1})
        continue
    end

    for icond = 1:2
        % create structures
        nt_cond = size(sims{isubj,icond}.rt,2);
        if isfirstpass
            idx_cond_sim = [zeros(1,nt_cond) ones(1,nt_cond)];
            resps_sim    = nan(nsim,nt_cond*2,nsubj);
            isfirstpass = false;
        end
        resps_sim(:,(1:nt_cond)+(nt_cond)*(icond-1),isubj) = sims{isubj,icond}.rt;
    end
end

out = struct;
out.resps   = resps_sim;
out.origdt  = sim_out.date; % date of sim generation
out.procdt  = datetime;     % date of processing 
out.nsim    = nsim;
out.pars    = sim_out.pars;

savedir = sprintf('../../processed/%s',samplename);
savename = sprintf('processed_%s',filename);
save(sprintf('%s/%s',savedir,savename),'out')