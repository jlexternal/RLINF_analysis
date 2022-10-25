% aggregate_output
clear all

samplename   = 'pilot07';
fittype      = 'recovery';
outputkernel = 'out_fit_noisyINF_recov_s'; % name kernel of output
momenttype   = 'xavg'; % xavg or xmap
npar = 3; % this might change depending on the type of model used
nsubj = 23;

load(sprintf('../constants/constants_rlinf_%s.mat',samplename),'ncnd'); % load constants
load(sprintf('../processed/%s/preprocessed_data_%s.mat',samplename,samplename),'idx_blmn');

pars = nan(nsubj,npar,ncnd);
out_vbmc = cell(nsubj,2);

for isubj = 1:nsubj
    % skip excluded subject indices
    if isnan(idx_blmn(isubj,1))
        continue
    end
    
    load(sprintf('./fitting/res/%s/%s%03d.mat',samplename,outputkernel,isubj));

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

% information about the data generated
out = struct;
out.samplename = samplename;
out.pars = pars;
out.outputkernel = outputkernel;
out.momenttype = momenttype;

% Write a description of the data file produced here
out.description = 'Fitted parameters on pilot07 with noisyINF model'; 

savename = 'pars_fit_rec_noisyINF.mat';
save(sprintf('./fitting/out/%s/%s',samplename,savename),'out');



