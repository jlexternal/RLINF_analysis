function [out] = fit_noisyINF_rlinf(cfg)

% check parameters
if ~isfield(cfg,'nstype')
    error('Missing learning noise type!');
elseif ~ismember('weber','white')
    error('Invalid learning noise type! (%s)',cfg.nstype)
end
if ~isfield(cfg,'fitalgo')
    cfg.fitalgo = '';
elseif ~ismember(cfg.fitalgo,{'bads','vbmc'})
    error('Invalid fitting algorithm! (%s)',cfg.fitalgo);
end
if ~isfield(cfg,'noprior')
    error('Missing priors flag!');
end
if ~isfield(cfg,'nsmp')
    error('Number of particle filter samples not provided!');
end
if ~isfield(cfg,'nres')
    cfg.nres = 1e2;
end
if ~isfield(cfg,'nrun')
    cfg.nrun = 1;
end
if ~isfield(cfg,'verbose')
    cfg.verbose = 0;
end
if ~isfield(cfg,'condstr')
    error('Missing condition type!')
elseif ~ismember(cfg.condstr,{'bandit','fairy'})
    error('Invalid condition type (%s)',cfg.condstr)
end

% get experiment data
if ~isfield(cfg,'cfg')
    % Note: experiment data should be trimmed/reshaped outside the fitting code!
    trl     = cfg.trl;  % trial number within block
    resp    = cfg.resp; % responses
    rt1     = cfg.rt1;   % reward/category strength of option 1
else
    trl     = cfg.cfg.trl;
    resp    = cfg.cfg.resp;
    rt1     = cfg.cfg.rt1;
end

% get total number of trials
ntrl = numel(trl);

% pass through fitting parameters
condstr = cfg.condstr; % condition type (bandit or fairy)
nstype  = cfg.nstype;  % noise type (weber or white)
fitalgo = cfg.fitalgo; % fitting algorithm (bads or vbmc)
noprior = cfg.noprior; % ignore priors?
nsmp    = cfg.nsmp;    % number of samples used by particle filter
nres    = cfg.nres;    % number of bootstrap/validation resamples
nrun    = cfg.nrun;    % number of random starting points (fitalgo = bads)
verbose = cfg.verbose; % fitting display level

% set internal parameters
epsi = 1e-6;
% set fixed statistics
m0     = 0.0000; % initial prior belief
llrmax = 1.8406; % llr of mostextreme stimulus
%!!! create a function that shows me how you've derived llrmax, so that we can be 100% sure of it!
% convert integer feedback values into llr values 
llr1 = conv_fb2llr(rt1,llrmax);


% define model parameters
np = 3; % number of parameters planned
pnam = cell(1,np); % name
pmin = nan(1,np);  % minimum value
pmax = nan(1,np);  % maximum value
pini = nan(1,np);  % initial value
pplb = nan(1,np);  % plausible lower bound
ppub = nan(1,np);  % plausible upper bound
pfun = cell(1,np);  % prior function (empty = uniform)
% 1/ hazard rate ~ uniform(0,1)
pnam{1} = 'h';
pmin(1) = 0.001;
pmax(1) = 0.999;
pini(1) = 0.5;
pplb(1) = 0.1;
ppub(1) = 0.9;
% 2/ inference noise
pnam{2} = 'sigma';
switch nstype
    case 'weber' % Weber noise ~ exp(1)
        pmin(2) = 0;
        pmax(2) = 30;
        pini(2) = 1;
        pplb(2) = 0.1;
        ppub(2) = 10;
        pfun{2} = @(x)exppdf(x,1);
    case 'white' % white noise ~ exp(0.1)
        pmin(2) = 0;
        pmax(2) = 10;
        pini(2) = 0.1;
        pplb(2) = 0.01;
        ppub(2) = 3;
        pfun{2} = @(x)exppdf(x,0.1);
end
% 3/ choice temperature ~ exp(0.1)
pnam{3} = 'tau';
pmin(3) = 1e-12;
pmax(3) = 10;
pini(3) = 1;
pplb(3) = 0.1;
ppub(3) = 5;
pfun{3} = @(x)exppdf(x,1);

ipar_noise = find(strcmpi(pnam,'sigma'));

if noprior
    % ignore priors
    pfun = cell(1,np);
end

% set number of parameters
npar = numel(pnam);

% apply user-defined initialization values
if isfield(cfg,'pini')
    for i = 1:npar
        if isfield(cfg.pini,pnam{i}) && ~isnan(cfg.pini.(pnam{i}))
            pini(i) = cfg.pini.(pnam{i});
            % clamp initialization value within plausible bounds
            pini(i) = min(max(pini(i),pplb(i)+1e-6),ppub(i)-1e-6);
        end
    end
end

% define fixed parameters
pfix = cell(1,npar);
for i = 1:npar
    if isfield(cfg,pnam{i}) && ~isempty(cfg.(pnam{i}))
        pfix{i} = min(max(cfg.(pnam{i}),pmin(i)),pmax(i));
    end
end

% define free parameters
ifit = cell(1,npar);
pfit_ini = [];
pfit_min = [];
pfit_max = [];
pfit_plb = [];
pfit_pub = [];
n = 1;
for i = 1:npar
    if isempty(pfix{i}) % free parameter
        ifit{i} = n;
        pfit_ini = cat(2,pfit_ini,pini(i));
        pfit_min = cat(2,pfit_min,pmin(i));
        pfit_max = cat(2,pfit_max,pmax(i));
        pfit_plb = cat(2,pfit_plb,pplb(i));
        pfit_pub = cat(2,pfit_pub,ppub(i));
        n = n+1;
    end
end
% set number of fitted parameters
nfit = length(pfit_ini);

% determine whether agent is noisy
isnoisy = isempty(pfix{ipar_noise}) || pfix{ipar_noise} > 0;

if nfit > 0
    switch fitalgo
        case 'bads'
            % fit model using Bayesian Adaptive Direct Search
            if ~exist('bads','file')
                error('BADS missing from path!');
            end
            if ~noprior
                % do not use priors
                warning('Disabling the use of priors when using BADS.');
                noprior = true;
            end
            % configure BADS
            options = bads('defaults');
            options.UncertaintyHandling = isnoisy; % noisy objective function
            options.NoiseFinalSamples = nres; % number of samples
            switch verbose % display level
                case 0, options.Display = 'none';
                case 1, options.Display = 'final';
                case 2, options.Display = 'iter';
            end
            % fit model using multiple random starting points
            fval   = nan(1,nrun);
            xhat   = cell(1,nrun);
            output = cell(1,nrun);
            for irun = 1:nrun
                done = false;
                while ~done
                    % set random starting point
                    n = 1;
                    for i = 1:npar
                        if isempty(pfix{i}) % free parameter
                            % sample starting point uniformly between plausible bounds
                            pfit_ini(n) = unifrnd(pplb(i),ppub(i));
                            n = n+1;
                        end
                    end

                    % fit model using BADS
                    % --------------------
                    [xhat{irun},fval(irun),exitflag,output{irun}] = ...
                        bads(@(x)getnl(x), ...
                        pfit_ini,pfit_min,pfit_max,pfit_plb,pfit_pub,[],options);
                    % --------------------

                    if exitflag > 0
                        done = true;
                    end
                end
            end
            % find best fit among random starting points
            [fval,irun] = min(fval);
            xhat   = xhat{irun};
            output = output{irun};
            % get best-fitting values
            phat = getpval(xhat);
            % create output structure with best-fitting values
            out = cell2struct(phat(:),pnam(:));
            % store fitting information
            out.fitalgo = fitalgo; % fitting algorithm
            out.nsmp    = nsmp;    % number of samples used by particle filter
            out.nres    = nres;    % number of validation resamples
            out.nrun    = nrun;    % number of random starting points
            out.ntrl    = ntrl;    % number of trials
            out.nfit    = nfit;    % number of fitted parameters
            % get maximum log-likelihood
            out.ll = -output.fval; % estimated log-likelihood
            out.ll_sd = output.fsd; % estimated s.d. of log-likelihood
            % get complexity-penalized fitting metrics
            out.aic = -2*out.ll+2*nfit+2*nfit*(nfit+1)/(ntrl-nfit+1); % AIC
            out.bic = -2*out.ll+nfit*log(ntrl); % BIC
            % get parameter values
            out.xnam = pnam(cellfun(@isempty,pfix));
            out.xhat = xhat;
            % store additional output from BADS
            out.output = output;
            
        case 'vbmc'
            % fit model using Variational Bayesian Monte Carlo
            if ~exist('vbmc','file')
                error('VBMC missing from path!');
            end
            % configure VBMC
            options = vbmc('defaults');
            options.MaxIter = 300; % maximum number of iterations
            options.MaxFunEvals = 500; % maximum number of function evaluations
            options.SpecifyTargetNoise = isnoisy; % noisy log-posterior function
            switch verbose % display level
                case 0, options.Display = 'none';
                case 1, options.Display = 'final';
                case 2, options.Display = 'iter';
            end
            
            % fit model using VBMC
            % --------------------
            [vp,elbo,~,exitflag,output] = vbmc(@(x)getlp(x), ...
                pfit_ini,pfit_min,pfit_max,pfit_plb,pfit_pub,options);
            % --------------------

            % generate 10^6 samples from the variational posterior
            xsmp = vbmc_rnd(vp,1e6);
            % get sample statistics
            xmap = vbmc_mode(vp);   % posterior mode
            xavg = mean(xsmp,1);    % posterior mean
            xstd = std(xsmp,[],1);  % posterior s.d.
            xcov = cov(xsmp);       % posterior covariance matrix
            xmed = median(xsmp,1);  % posterior medians
            xqrt = quantile(xsmp,[0.25,0.75],1); % posterior 1st and 3rd quartiles
            % get full parameter set with best-fitting values
            phat_map = getpval(xmap); % posterior mode
            phat_avg = getpval(xavg); % posterior mean
            % use posterior mode as default
            phat = phat_map;
            
            % create output structure
            out = cell2struct(phat(:),pnam(:));
            % create substructure with posterior mode
            out.pmap = cell2struct(phat_map(:),pnam(:));
            % create substructure with posterior mean
            out.pavg = cell2struct(phat_avg(:),pnam(:));
            % store fitting information
            out.fitalgo = fitalgo; % fitting algorithm
            out.nsmp    = nsmp;    % number of samples used by particle filter
            out.nres    = nres;    % number of bootstrap resamples
            out.ntrl    = ntrl;    % number of trials
            out.nfit    = nfit;    % number of fitted parameters
            % store variational posterior solution
            out.vp = vp;
            % get ELBO (expected lower bound on log-marginal likelihood)
            out.elbo = elbo; % estimate
            out.elbo_sd = output.elbo_sd; % standard deviation
            % get maximum log-posterior and maximum log-likelihood
            out.lp = getlp(xmap); % log-posterior
            out.ll = getll(phat_map{:}); % log-likelihood
            % get parameter values
            out.xnam = pnam(cellfun(@isempty,pfix)); % fitted parameters
            out.xmap = xmap; % posterior mode
            out.xavg = xavg; % posterior mean
            out.xstd = xstd; % posterior s.d.
            out.xcov = xcov; % posterior covariance matrix
            out.xmed = xmed; % posterior median
            out.xqrt = xqrt; % posterior 1st and 3rd quartiles
            % store extra VBMC output
            out.output = output;
            
        otherwise
            error('Undefined fitting algorithm!');   
    end
    
else
    % use fixed parameter values
    phat = getpval([]);
    
    % create output structure
    out = cell2struct(phat(:),pnam(:));
    
    % run particle filter
    [pt_hat,lt_hat,ut_hat,st_hat,wt_hat] = getp(phat{:});
    
    % average trajectories
    pt_avg = mean(pt_hat,2);
    mt_avg = sum(bsxfun(@times,mt_hat,reshape(wt_hat,[ntrl,1,nsmp])),3);
    ut_avg = sum(bsxfun(@times,ut_hat,reshape(wt_hat,[ntrl,1,nsmp])),3);
    
    % store averaged trajectories
    out.pt = pt_avg; % response probabilities
    out.mt = mt_avg; % filtered posterior means
    out.vt = vt_hat; % posterior variances
    
    % DEBUG: needs rework
    % identify greedy responses
    [~,ru] = max(ut_avg,[],2); % based on exact updates
    [~,rf] = max(mt_avg,[],2); % based on noisy updates
    % compute fractions of non-greedy responses
    pngu = mean(resp ~= ru); % based on exact updates
    pngf = mean(resp ~= rf); % based on noisy updates
    pngd = mean(rf ~= ru); % due to noise
    pngx = mean(resp(resp ~= ru) == rf(resp ~= ru)); % explained by noise

    % store fractions of non-greedy responses
    out.pngu = pngu; % based on exact updates
    out.pngf = pngf; % based on noisy updates
    out.pngd = pngd; % due to noise
    out.pngx = pngx; % explained by noise
    
end

% store configuration structure
out.cfg = cfg;

    %-------------------------------------
    function [pval] = getpval(p)
        % get parameter values
        pval = cell(1,npar);
        for k = 1:npar
            if isempty(pfix{k}) % free parameter
                pval{k} = p(ifit{k});
            else % fixed parameter
                pval{k} = pfix{k};
            end
        end
    end
    %-------------------------------------
    function [nl] = getnl(p)
        pval = getpval(p);      % get parameter values
        nl = -getll(pval{:});   % get negative log-likelihood
    end
    %-------------------------------------
    function [lp,lp_sd] = getlp(p)
        pval = getpval(p); % get parameter values
        % get log-prior
        l0 = 0;
        for k = 1:npar
            if isempty(pfix{k}) % free parameter
                if isempty(pfun{k}) % use uniform prior
                    l0 = l0+log(unifpdf(pval{k},pmin(k),pmax(k)));
                else % use specified prior
                    l0 = l0+log(pfun{k}(pval{k}));
                end
            end
        end
        [ll,ll_sd] = getll(pval{:}); % get log-likelihood
        % get log-posterior
        lp = ll+l0; % estimate
        lp_sd = ll_sd; % bootstrap s.d.
    end
    %-------------------------------------
    function [ll,ll_sd] = getll(varargin)
        % compute response probability
        p = getp(varargin{:});
        if nargout > 1
            % compute log-likelihood s.d.
            lres = nan(nres,1);
            for ires = 1:nres
                jres = randsample(nsmp,nsmp,true);
                pres = mean(p(:,jres),2);
                pres = epsi+(1-epsi*2)*pres;
                lres(ires) = ...
                    sum(log(pres(resp == 1)))+ ...
                    sum(log(1-pres(resp == 2)));
            end
            ll_sd = max(std(lres),1e-6);
        end
        % compute log-likelihood
        p = mean(p,2);
        p = epsi+(1-epsi*2)*p;
        ll = ...
            sum(log(p(resp == 1)))+ ...
            sum(log(1-p(resp == 2)));
    end
    %-------------------------------------
    function [pt,lt,ut,st,wt] = getp(h,sigma,tau)

        % initialize output variables
        pt = nan(ntrl,nsmp); % probability of blue/moon (option 1)
        lt = nan(ntrl,nsmp); % option 1 relative belief
        ut = nan(ntrl,nsmp); % unfiltered option 1 relative belief
        et = nan(ntrl,nsmp); % prediction errors
        st = nan(ntrl,nsmp); % filtering noise
        wt = nan(ntrl,nsmp); % filtering weights
        
        % run filter
        for itrl = 1:ntrl
            if trl(itrl) == 1
                % initialize posterior means and variance
                lt(itrl,:) = m0;
                % respond randomly
                if strcmpi(condstr,'bandit')
                    pt(itrl,:) = 0.5;
                end
                % initialize filtering noise and weights
                st(itrl,:,:) = 0;
                wt(itrl,:) = 1/nsmp;
                continue
            end
            
            % condition posterior means using weighted bootstrapping
            lt(itrl,:)  = lt(itrl-1,randsample(nsmp,nsmp,true,wt(itrl-1,:)));

            % prediction error sources
            if strcmpi(condstr,'bandit')
                if isnan(llr1(itrl-1))
                    error('NaN on llr1 at itrl-1 = %d',itrl-1);
                end
                %!!! let's not use this prediction error term
                %!!! if anything, we should use the KL-divergence between prior and posterior beliefs
                %!!! like KLdiv(lt(itrl,:),psi(lt(itrl,:),h) + llr1(itrl-1))
                et(itrl,:) = llr1(itrl-1)-lt(itrl-1,:); % "prediction error"
                %!!! below psi(lt(itrl-1),h) is wrong because you are using the non-filtered posterior belief
                %!!! this should be ut(itrl,:) = psi(lt(itrl,:),h) + llr1(itrl-1);
                ut(itrl,:) = psi(lt(itrl-1),h) + llr1(itrl-1);
            else
                if isnan(llr1(itrl)) %debug
                    error('NaN on rt1 at itrl = %d',itrl);
                end
                % Fairy: stimulus(t) - *learning(t)* - choice
                %!!! same as above: let's not use this prediction error term
                %!!! this could be in this case KLdiv(lt(itrl,:),psi(lt(itrl,:),h) + llr1(itrl))
                et(itrl,:) = llr1(itrl)-lt(itrl,:);
                %!!! same as above: below psi(lt(itrl-1),h) is wrong because you are using the non-filtered posterior belief
                %!!! this should be ut(itrl,:) = psi(lt(itrl,:),h) + llr1(itrl);
                ut(itrl,:) = psi(lt(itrl-1),h) + llr1(itrl);
            end
            
            % calculate std of error 
            if strcmp(nstype,'weber')
                %!!! let's not use this one, unless you use the KLdiv() definition of the prediction error
                st(itrl,:) = sigma*abs(et(itrl,:));
            else
                %!!! this is fine!
                st(itrl,:) = sigma;
            end

            % account for filtering noise
            if isnoisy
                %!!! this is fine!
                lt(itrl,:) = normrnd(ut(itrl,:),st(itrl,:));
            end
            
            % apply policy to get response probabilities of option 1
            %!!! it would actually be good to use a sigma-like definition of choice stochasticity
            %!!! pt(itrl,:) = normcdf(lt(itrl,:),0,sigma_sel);
            %!!! where sigma_sel is the amount of response selection noise
            %!!! you would then rename sigma as sigma_inf for inference noise
            x = reshape(lt(itrl,:),[1,nsmp]);
            pt(itrl,:) = 1./(1+exp(-x/tau));

            % compute filtering weights
            if resp(itrl) == 1
                wt(itrl,:) = pt(itrl,:);
            else
                wt(itrl,:) = 1-pt(itrl,:);
            end
            if nnz(wt(itrl,:)) == 0
                wt(itrl,:) = 1/nsmp;
            else
                wt(itrl,:) = wt(itrl,:)/sum(wt(itrl,:));
            end
        end
    end
    %-------------------------------------
    function llr = conv_fb2llr(fb,llrmax)
        % Takes the feedback value (task variable) relative to the true state and converts it
        % to the LLR value (evidence) given the generative FNR (false negative rate)
        lvl_list = (-49:+49)/49; % list of color levels
        llr_list = lvl_list*llrmax; % list of logLRs
        llr_list(end+1) = nan;
        ifinal = numel(llr_list);
        fb(isnan(fb)) = ifinal;
        llr = llr_list(fb);
    end
    %-------------------------------------
    function lprev = psi(l,h)
        % model function that transforms prior based on hazard rate (Glaze et al. 2015)
        lprev = l + log((1-h)/h + exp(-l)) - log((1-h)/h + exp(l));
    end
    %-------------------------------------
end