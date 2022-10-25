% compare_fits
clear all

samplename  = 'pilot07';
kernelarr   = {'out_fit_noisyINF_ALL',...
               'out_fit_noisyINF_zeroinf_ALL',...
               'out_fit_noisyINF_zerosel_ALL'};
momenttype   = 'xavg'; % xavg or xmap
nmod = numel(kernelarr);

nsubj = 23;

load(sprintf('../../constants/constants_rlinf_%s.mat',samplename),'ncnd'); % load constants
load(sprintf('../../processed/%s/preprocessed_data_%s.mat',samplename,samplename),'idx_blmn');

elbos = nan(nsubj,ncnd,nmod);

for imod = 1:nmod
    load(sprintf('./sample_in/%s/%s.mat',samplename,kernelarr{imod}));
    for isubj = 1:nsubj
        if isempty(out_vbmc{isubj,1})
            continue
        end
        for icond = 1:ncnd
            elbos(isubj,icond,imod) = out_vbmc{isubj,icond}.elbo;
        end
    end
end

