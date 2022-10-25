function [out] = sim_noisyKF_rlinf(cfg)

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
if ~all(isfield(cfg,{'alpha','zeta','tau','nstype'}))
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
nstype  = cfg.nstype;
alpha   = cfg.alpha;
zeta    = cfg.zeta;
tau     = cfg.tau;
nsim    = cfg.nsim;

isnoisy = zeta ~= 0;

% set fixed statistics
m0 = 0.0000; % prior mean
v0 = 0.0214; % prior variance
vs = 0.0643; % sampling variance

% define reparameterization functions:
%   * alpha = fa(vd/vs)
%   * vd/vs = fv(alpha)
fa = @(v)1./(1+exp(+0.4486-log2(v)*0.6282)).^0.5057;
av = 0.001:0.001:0.999;
vv = arrayfun(@(a)fzero(@(v)fa(v)-a,2.^[-30,+30]),av);
fv = @(a)interp1(av,vv,a,'pchip');

% clip parameter values for numerical stability
alpha = min(max(alpha,0.001),0.999);
tau   = max(tau,1e-12);

% set drift variance 
vd = fv(alpha)*vs;

% tracked variables
rt = nan(nsim,nt); % responses
pt = nan(nsim,nt); % probability of blue/moon (option 1)
mt = nan(nsim,nt); % option 1 relative values
ut = nan(nsim,nt); % unfiltered option 1 relative values
et = nan(nsim,nt); % error terms
st = nan(nsim,nt); % filtering noise
vt = nan(1,nt);  % posterior variances

ib = 0; % block counter
for it = 1:nt
    % initialize for beginning of block
    if trl(it) == 1
        ib = ib + 1;
        mt(:,it) = m0;
        vt(:,it) = v0;
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
        et(:,it)  = (r1(it-1)-.5)-mt(:,it-1); 
    else
        et(:,it)  = (r1(it)-.5)-mt(:,it-1);
    end
    % compute Kalman gain
    kgain = vt(it-1)/(vt(it-1)+vs);

    % update posterior means and variances
    ut(:,it) = mt(:,it-1) + kgain*et(:,it);
    if strcmp(nstype,'weber')
        st(:,it) = zeta*kgain*abs(et(:,it));
    else
        st(:,it) = zeta;
    end
    mt(:,it) = ut(:,it);
    if isnoisy
        mt(:,it) = noisrnd(ut(:,it),st(:,it));
    end
    vt(it)   = (1-kgain)*vt(it-1); % Kalman update on posterior variance
    vt(it)   = vt(it) + vd; % account for drift

    % choice step
    pt(:,it) = 1./(1+exp(-mt(:,it)/tau));
    rt(:,it) = (rand(nsim,1) > pt(:,it))+1;
end

out = [];
out.m0    = m0;    % prior mean
out.v0    = v0;    % prior variance
out.vs    = vs;    % sampling variance
out.vd    = vd;    % drift variance
out.alpha = alpha; % learning rate
out.zeta  = zeta;  % learning noise
out.tau   = tau;   % policy temperature

out.pt    = pt;    % response probabilities
out.rt    = rt;    % responses
out.ut    = ut;    % unfiltered posterior means
out.mt    = mt;    % posterior means
out.vt    = vt;    % posterior variances


end
% -------------------------------------------------------------
function [x] = noisrnd(m,s)
    % sample from noise distribution = truncated normal in [0,1]
    if numel(s) == 1
        s = s(ones(size(m)));
    end
    x = nan(size(m));
    f = true(size(m));
    while any(f(:))
        x(f) = normrnd(m(f),s(f));
        f(x >= -.5 & x <= .5) = false;
    end
end