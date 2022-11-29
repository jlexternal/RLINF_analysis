function [out] = sim_noisyKF_cfrule_rlinf(cfg)
% check input arguments
if isfield(cfg,'getInfo')
    disp('Description:');
    fprintf(['Simulation function for the 2 option noisy Kalman filter with' ...
             ' counterfactual rule toggle.\n\n'])
    disp('Input parameters (config structure):');
    disp(' r1:        Option 1 value data     [0 1]');
    disp(' trl:       Trial number            (+''ve integers)');
    disp(' nsim:      Number of simulations   (+''ve integer)');
    disp(' condstr:   Condition string        (''bandit''/''fairy'')');
    disp(' firstresp: First responses         {1,2}');
    fprintf('\n');
    disp('Model parameters:');
    disp('  name   description    range');
    disp(' nstype: Learning noise type. (''weber'' or ''white'')');
    disp(' cfrule: Counterfactual rule. (true/false)');
    disp(' alpha:  Learning rate.  [0 1]');
    disp(' delta:  Decay rate.     [0 1]');
    disp(' zeta:   Learning noise. [0 10)');
    disp(' tau:    Choice temp.    (0 5)');
    return
end
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
if ~all(isfield(cfg,{'alpha','zeta','tau','delta','nstype','cfrule'}))
    error('Missing model parameters!');
end
if ~isfield(cfg,'nsim')
    error('Missing number of simulations!');
end
rand1st = true; % random first response flag
if isfield(cfg,'firstresp')
    rand1st = false;
    resp1 = cfg.firstresp; % (vector) first trial responses for each block
end

%---------
% get experiment information
condstr = cfg.condstr;  % (string) condition type
r1      = cfg.r1;       % (vector) value of blue/moon
if size(r1,2) ~= 1
    r1 = r1';
end
trl     = cfg.trl;      % (vector) trial index
nt      = numel(r1);    % number of trials
r2      = 1-r1;         % (vector) value of NOT blue/moon
rew     = [r1 r2];      % 2-column matrix of reward values for both options

% localize input parameters
cfrule  = cfg.cfrule;
nstype  = cfg.nstype;
alpha   = cfg.alpha;
delta   = cfg.delta;
zeta    = cfg.zeta;
tau     = cfg.tau;
nsim    = cfg.nsim;

isnoisy = zeta ~= 0;

% set fixed statistics
m0 = 0.5000; % prior mean
v0 = 0.0214; % prior variance
vs = 0.0163; % sampling variance

% define reparameterization functions:
%   * alpha = fa(vd/vs)
%   * vd/vs = fv(alpha)
fa = @(v)1./(1+exp(+0.4486-log2(v)*0.6282)).^0.5057;
av = 0.001:0.001:0.999;
vv = arrayfun(@(a)fzero(@(v)fa(v)-a,2.^[-30,+30]),av);
fv = @(a)interp1(av,vv,a,'pchip');

% clip parameter values to avoid numerical instability
alpha = min(max(alpha,0.001),0.999);
tau   = max(tau,1e-12);

% get drift variance
vd = fv(alpha)*vs;

% tracked variables
rt = nan(nsim,nt);   % responses
pt = nan(nsim,nt);   % probability of blue/moon (option 1)
mt = nan(nsim,nt,2); % filtered posterior means 
ut = nan(nsim,nt,2); % unfiltered posterior means 
et = nan(nsim,nt,2); % error terms
st = nan(nsim,nt,2); % filtering noise
vt = nan(nsim,nt,2);  % posterior variances

r_ch = nan(nsim,nt); % chosen rewards
r_un = nan(nsim,nt); % unchosen rewards

ib = 0; % block counter
for it = 1:nt
    % initialize block (first trial)
    if trl(it) == 1
        ib = ib + 1;
        % initialize posterior means and variances
        mt(:,it,:) = m0;
        vt(:,it,:) = v0;
        pt(:,it)  = 0.5; % flat prior
        st(:,it)  = 0;

        % initialize bandit first response
        if strcmpi(condstr,'bandit') % bandit task
            if rand1st
                rt(:,it) = round(rand(nsim,1)) + 1;
            else
                rt(:,it) = resp1(ib)*ones(nsim,1);
            end
            % get chosen and unchosen rewards
            r_ch(:,it) = rew(sub2ind(size(rew),it(ones(nsim,1)),rt(:,it)));
            r_un(:,it) = rew(sub2ind(size(rew),it(ones(nsim,1)),3-rt(:,it)));
        else % fairy task (first response is nothing)
            rt(:,it) = nan(nsim,1);
        end
        continue
    end

    % compute Kalman gain
    kgain = vt(:,it-1,:)./(vt(:,it-1,:)+vs);

    % initialize fairy first response
    if strcmpi(condstr,'fairy') & trl(it) == 2
        et(:,it,1) = r1(it)-mt(:,it-1,1);
        et(:,it,2) = (1-r1(it))-mt(:,it-1,2);
        ut(:,it,:) = mt(:,it-1,:)+kgain(:,1,:).*et(:,it,:);
        mt(:,it,:) = ut(:,it,:);
        vt(:,it,:) = (1-kgain(:,1,:)).*vt(:,it-1,:);
        if strcmpi(nstype,'weber')
            st(:,it,:) = zeta*kgain(:,1,:).*abs(et(:,it,:));
        else
            st(:,it,:) = zeta;
        end
        if isnoisy
            mt(:,it,:) = noisrnd_cf(ut(:,it,:),st(:,it,:));
        end
        pt(:,it) = 1./(1+exp(-(mt(:,it,1)-mt(:,it,2))/tau));

        if ~rand1st 
            rt(:,it) = resp1(ib)*ones(nsim,1);
        else
            rt(:,it) = (rand(nsim,1) > pt(:,it))+1;
        end
        continue
    end

    % -------------------------------------------------------------- %
    % prediction error for both options where chosen option = 1
    i1 = rt(:,it-1) == 1;
    et(i1,it,1)  = r1(it-1)-mt(i1,it-1,1); 
    et(i1,it,2)  = (1-r1(it-1))-mt(i1,it-1,2); 

    % update posterior mean and variance of chosen option
    ut(i1,it,1) = mt(i1,it-1,1)+kgain(i1,1,1).*et(i1,it,1);
    % corrupt tracked value with noise
    if strcmpi(nstype,'weber')
        st(i1,it,1) = zeta*kgain(i1,1,1).*abs(et(i1,it,1));
    else
        st(i1,it,1) = zeta;
    end
    mt(i1,it,1) = ut(i1,it,1);
    if isnoisy
        mt(i1,it,1) = noisrnd_cf(ut(i1,it,1),st(i1,it,1));
    end
    
    vt(i1,it,1) = (1-kgain(i1,1,1)).*vt(i1,it-1,1);
    
    if ~cfrule
        % decay posterior mean of unchosen option
        mt(i1,it,2) = mt(i1,it-1,2)+delta*(0.5-mt(i1,it-1,2));
        ut(i1,it,2) = mt(i1,it,2);
        vt(i1,it,2) = vt(i1,it-1,2);
    else
        % update posterior mean of unchosen option using counterfactual rule
        ut(i1,it,2) = mt(i1,it-1,2)+kgain(i1,1,2).*et(i1,it,2);
        if strcmpi(nstype,'weber')
            st(i1,it,2) = zeta*kgain(i1,1,2).*abs(et(i1,it,2));
        else
            st(i1,it,2) = zeta;
        end
        mt(i1,it,2) = ut(i1,it,2);
        if isnoisy
            mt(i1,it,2) = noisrnd_cf(ut(i1,it,2),st(i1,it,2));
        end
        vt(i1,it,2) = (1-kgain(i1,1,2)).*vt(i1,it-1,2);
    end

    % -------------------------------------------------------------- %
    % prediction error for both options where chosen option = 2
    i2 = rt(:,it-1) == 2 & ~isnan(rew(it-1,1));
    et(i2,it,2)  = (1-r1(it-1))-mt(i2,it-1,2); 
    et(i2,it,1)  = r1(it-1)-mt(i2,it-1,1);
    
    % update posterior mean and variance of chosen option
    ut(i2,it,2) = mt(i2,it-1,2)+kgain(i2,1,2).*et(i2,it,2);
    if strcmpi(nstype,'weber')
        st(i2,it,2) = zeta*kgain(i2,1,2).*abs(et(i2,it,2));
    else
        st(i2,it,2) = zeta;
    end
    mt(i2,it,2) = ut(i2,it,2);
    if isnoisy
        mt(i2,it,2) = noisrnd_cf(ut(i2,it,2),st(i2,it,2));
    end
    vt(i2,it,2) = (1-kgain(i2,1,2)).*vt(i2,it-1,2);
    
    if ~cfrule
        % decay posterior mean of unchosen option
        mt(i2,it,1) = mt(i2,it-1,1)+delta*(0.5-mt(i2,it-1,1));
        ut(i2,it,1) = mt(i2,it,1);
        vt(i2,it,1) = vt(i2,it-1,1);
    else
        % update posterior mean of unchosen option using counterfactual rule
        ut(i2,it,1) = mt(i2,it-1,1)+kgain(i2,1,1).*et(i2,it,1);
        if strcmpi(nstype,'weber')
            st(i2,it,1) = zeta*kgain(i2,1,1).*abs(et(i2,it,1));
        else
            st(i2,it,1) = zeta;
        end
        mt(i2,it,1) = ut(i2,it,1);
        if isnoisy
            mt(i2,it,1) = noisrnd_cf(ut(i2,it,1),st(i2,it,1));
        end
        vt(i2,it,1) = (1-kgain(i2,1,1)).*vt(i2,it-1,1);
    end
    
    % account for drift
    vt(:,it,:) = vt(:,it,:)+vd;
    
    % compute log-likelihood and log-variance ratios
    llr(:,it) = 2*(mt(:,it,1)-mt(:,it,2))./sqrt(vt(:,it,1)+vt(:,it,2));
    lvr(:,it) = log(vt(:,it,1))-log(vt(:,it,2));
    
    % apply policy to choose response
    pt(:,it) = 1./(1+exp(-(mt(:,it,1)-mt(:,it,2))/tau)); 
    % get response
    rt(:,it) = (rand(nsim,1) > pt(:,it))+1;
    % get response properties
    rgrd(:,it) = rt(:,it) == 1+(llr(:,it) < 0); % greediness
    rcur(:,it) = rt(:,it) == 1+(lvr(:,it) < 0); % curiosity
    % get chosen and unchosen rewards
    r_ch(:,it) = rew(sub2ind(size(rew),it(ones(nsim,1)),rt(:,it)));
    r_un(:,it) = rew(sub2ind(size(rew),it(ones(nsim,1)),3-rt(:,it)));
    
end

% remove respones in response matrix where there was no trial
resp(isnan(rt(:,:,1))) = nan;

% create output structure
out       = [];
out.m0    = m0;    % prior mean
out.v0    = v0;    % prior variance
out.vs    = vs;    % sampling variance
out.vd    = vd;    % drift variance
out.alpha = alpha; % learning rate
out.delta = delta; % learning decay
out.zeta  = zeta;  % learning noise
out.tau   = tau;   % policy temperature
out.cfrule = cfrule; % counterfactual rule
out.nstype = nstype; % noise type

out.rt    = rt;    % responses
out.pt    = pt;    % response probabilities
out.ut    = ut;    % unfiltered posterior means
out.mt    = mt;    % posterior means
out.vt    = vt;    % posterior variances

out.rgrd  = rgrd;  % response greediness
out.rcur  = rcur;  % response curiosity
out.r_ch  = r_ch;  % obtained rewards
out.r_un  = r_un;  % foregone rewards
out.llr   = llr;   % log-likelihood ratio
out.lvr   = lvr;   % log-variance ratio

end