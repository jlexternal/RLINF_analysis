% compare_KF_lesions
clear all

datatypes = {'full','nozeta','notau'};
lesiontypes = {'none','zeta','tau'};
simdattypes = {'','nozeta_','notau_'};

icond = 1;

elbos = nan(100,3,3); % sets, model(lesions), data source
idx_win = zeros(100,3,3);

pars = zeros(100,3,3,3); % sets, parameters, model, data source

pars_gen = nan(100,3,3); % sets, parameters, data source

for iset = 1:100
    for idat = 1:3
        % load generative parameters for simulation data
        filename = sprintf('sim_dat_4confusion_%sKF.mat',simdattypes{idat});
        load(filename);
        pars_gen(iset,:,idat) = [sim.dat{iset,icond}.alpha sim.dat{iset,icond}.zeta sim.dat{iset,icond}.tau];


        for iles = 1:3
            filename = sprintf('out_fit_noisyKF_conf_%sData_%sLesionedModel_s%03d.mat',datatypes{idat},lesiontypes{iles},iset);
            load(sprintf('./fits/%s',filename));

            elbos(iset,iles,idat) = out_fit.out_vbmc{icond}.elbo;
            switch iles
                case 1
                    parset = 1:3;
                case 2
                    parset = [1 3];
                case 3
                    parset = 1:2;
            end
            pars(iset,parset,iles,idat) = out_fit.out_vbmc{icond}.xavg;

        end
        [~,i] = max(elbos(iset,:,idat));
        idx_win(iset,i,idat) = 1;
    end
end

%%
modelfreq = squeeze(sum(idx_win))'/100
figure(1)
imagesc(modelfreq)
hold on
colorbar

%%
isource = 1;
imodel = 1;
[r,p] = corr(pars_gen(:,:,isource),pars(:,:,imodel,isource))
figure(2)
clf
for i = 1:4
    subplot(1,4,i);
    if i < 4
        x = pars_gen(:,i,isource);
        y = pars(:,i,imodel,isource);
        scatter(x,y);
        hold on
        plot([0 max([x,y])],[0 max([x,y])])
        xlim([0,max(y)]);
        ylim([0 max(y)]);
    else
        imagesc(r)
    end

end