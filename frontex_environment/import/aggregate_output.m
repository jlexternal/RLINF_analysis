% aggregate_output
clear all

% INPUT ---------------------------------
samplename   = 'sample1';
fittype      = 'fit'; % 'fit' or 'recovery'
outputkernel = 'out_fit_noisyKF_s'; % name kernel of output
momenttype   = 'xavg'; % xavg or xmap
npar = 3; % this might change depending on the type of model used
nsubj = 113;
% --------------------------------------

load(sprintf('../../constants/constants_rlinf_%s.mat',samplename),'ncnd'); % load constants
load(sprintf('../../processed/%s/preprocessed_data_%s.mat',samplename,samplename),'idx_blmn');

pars = nan(nsubj,npar,ncnd);
out_vbmc = cell(nsubj,2);

for isubj = 1:nsubj
    % skip excluded subject indices
    if isnan(idx_blmn(isubj,1))
        continue
    end
    
    load(sprintf('./sample_in/%s/%s%03d.mat',samplename,outputkernel,isubj));

    switch fittype
        case 'fit'
            out_vbmc(isubj,:) = out_fit.out_vbmc(isubj,:);
        case 'recovery'
            out_vbmc(isubj,:) = out_rec.out_vbmc(isubj,:);
    end

    for icond = 1:2
        pars(isubj,:,icond) = out_vbmc{isubj,icond}.(momenttype);
    end
end

% save the aggregated vbmc output structures with the rest of the
% individual fits
save(sprintf('./sample_in/%s/%s_ALL.mat',samplename,outputkernel(1:end-2)),'out_vbmc');

% information about the data generated
out = struct;
out.samplename = samplename;
out.pars = pars;
out.outputkernel = outputkernel;
out.momenttype = momenttype;

% Write a description of the data file produced here
out.description = 'Recovered parameters on sample1 (subj_max 113) with noisyKF model'; 

savename = 'pars_fit_noisyKF_ALL.mat';
save(sprintf('./sample_out/%s/%s',samplename,savename),'out');



