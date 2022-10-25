function [out] = sim_noisyINF_rlinf(cfg)

% check input arguments
if ~isfield(cfg,'condstr')
    error('Missing condition string! (condstr)');
    if ~ismember(cfg.condstr,{'bandit','fairy'})
        error('Unexpected value for condstr (''%s'')!',cfg.condstr);
    end
end
if ~isfield(cfg,'r1')
    error('Missing option value data! (r1)!');
end
if ~isfield(cfg,'trl')
    error('Missing trial number array! (trl)');
end
if ~all(isfield(cfg,{'h','sigma_inf','sigma_sel'}))
    error('Missing model parameters!');
end
if ~isfield(cfg,'nsim')
    error('Missing number of simulations!');
end
rand1st = true; % random first response flag
if isfield(cfg,'resp1')
    rand1st = false;
    resp1 = cfg.firstresp; % (vector) first trial responses for each block
end

% get experiment information
condstr = cfg.condstr;   % (string) condition type
r1      = cfg.r1;        % (vector) value of blue/moon
trl     = cfg.trl;       % (vector) trial index
nt      = numel(r1);     % number of trials

% localize input parameters
h         = cfg.h;
sigma_inf = cfg.sigma_inf;
sigma_sel = cfg.sigma_sel;
nsim      = cfg.nsim;

% set fixed statistics
l0 = 0.0000; % prior mean
llrmax = 1.8406; % llr of mostextreme stimulus

% convert feedback values to LLRs
llr1    = conv_fb2llr(r1,llrmax);

% clip parameter values for numerical stability
h         = min(max(h,0.001),0.999);
sigma_sel = max(sigma_sel,1e-12);

% tracked variables
rt = nan(nsim,nt); % responses
ut = nan(nsim,nt); % unfiltered option 1 relative belief
lt = nan(nsim,nt); % option 1 relative belief
st = nan(nsim,nt); % filtering noise
pt = nan(nsim,nt); % probability of blue/moon (option 1)

ib = 0; % block counter
for it = 1:nt
    % initialize for beginning of block
    if trl(it) == 1
        ib = ib + 1;
        lt(:,it) = l0;
        pt(:,it) = 0.5; % flat prior
        st(:,it) = 0;
        if strcmpi(condstr,'bandit') % bandit task
            if rand1st
                rt(:,it) = round(rand(nsim,1)) + 1;
            else
                rt(:,it) = repmat(resp1(ib)*ones(nsim,1));
            end
        else % fairy task (first response is nothing)
            rt(:,it) = nan(nsim,1);
        end
        continue
    end
    
    % prediction error 
    if strcmpi(condstr,'bandit') % bandit task
        ut(:,it) = psi(lt(:,it-1),h) + llr1(it-1);
    else
        ut(:,it) = psi(lt(:,it-1),h) + llr1(it);
    end
    
    st(:,it) = sigma_inf;
    % corrupt belief
    lt(:,it) = normrnd(ut(:,it),st(:,it));

    % choice step
    pt(:,it) = normcdf(lt(:,it),0,sigma_sel);
    rt(:,it) = (rand(nsim,1) > pt(:,it))+1;
end

out = [];
out.l0        = l0;    % prior mean
out.h         = h;     % hazard rate
out.sigma_inf = sigma_inf; % inference noise
out.sigma_sel = sigma_sel;   % policy temperature

out.pt    = pt;    % response probabilities
out.rt    = rt;    % responses
out.ut    = ut;    % unfiltered posterior belief

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