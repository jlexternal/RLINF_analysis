function simulate_kf(cfg)
% Description: Simulates the KF on some instance of the RLINF task
%
% Input:
%   cfg:
%   .task: arrays of task information
%   .nsim: number of simulations
%   .pars: parameters
%   .resp: subject responses 
%   .rnd1: random first choice
%   .hyperpars: hyperparameters if not default

% Check for configuration variables
if ~isfield(cfg,'task')
    error('Task data missing in config structure!');
end
% Number of simulation samples
if ~isfield(cfg,'nsim')
    nsim = 1e4;
else
    nsim = cfg.nsim;
end
% Noise scaling type
if ~isfield(cfg,'nstype')   % weber or white 
    nstype = cfg.nstype;
else
    nstype = 'weber';       % default Weber noise scaling
end
% Choice policy
if ~isfield(cfg,'chrule')   % softm, argm, thomp, ucb
    chrule = cfg.chrule;
else
    chrule = 'argmax';      % default Argmax choice policy
end
if ~isfield(cfg,'params')
    error('KF parameters missing from config structure!');
end
% random first choice
rnd1 = true;
if ~isfield(cfg,'rnd1')
    rnd1 = cfg.rnd1;
end

% extract task data
r1 = cfg.task.r1; % value option 1
r2 = cfg.task.r2; % value option 2
idx_trl = cfg.task.trl; % trial indices
idx_blk = cfg.task.blk; % block indices
isrl  = cfg.task.isrl;  % whether the trial belongs to RL task (or INF)

% Set hyperparameters/fixed statistics of KF (maybe unnecessary)
m0 = 0.5000; % prior mean
if isfield(cfg,'hyperpars')
    v0 = cfg.hyperpars.v0;
    vs = cfg.hyperpars.vs;
else
    v0 = 0.0214; % prior variance
    vs = 0.0163; % sampling variance 
end

% Instantiate model
if rnd1 
    fprintf('Assuming random first choice for each block...\n');
else
    fprintf('Matching simulation responses to subject''s first choice...\n');
end

nt = numel(r1);     % number of trials in task
nb = max(idx_blk);  % number of blocks in task

% default values for KF parameters
alpha   = .5;
zeta    = 0;
tau     = 0;

if isvar(cfg.params.alpha)
    alpha = cfg.params.alpha;  % learning rate / perceived volatility
end

if ~exists(cfg.params.zeta)
    zeta = cfg.params.zeta;    % learning noise
else
    zeta = 0;
end
tau = cfg.tau;      % choice temperature

if zeta == 0 & tau == 0
    fprintf('No noise added. Setting number of samples to 1..\n');
    nsim = 1;
end

% tracked values
kt = nan(nt,2,nsim); % Kalman gain value
mt = nan(nt,2,nsim); % posterior means
vt = nan(nt,2,nsim); % posterior variances
ut = nan(nt,2,nsim); % filtering means
st = nan(nt,2,nsim); % filtering noise
et = nan(nt,2,nsim); % prediction errors
pt = nan(nt,1,nsim); % tracked option choice probability
rt = nan(nt,1,nsim); % choices made

% Task structure
% RL task:
%   Choice(t-1) > Result(t) > Update(t) > Choice(t)
%
% Inference task:
%   Result(t) > Update(t) > Choice(t+1) > Result(t+1)

% DEBUG: need to sync the two types of tasks together

fprintf('hello world')

for it = 1:nt
    ib = idx_blk(it); % block index

    % initialize/reset KF
    if idx_trl(it) == 1 % first trial of block
        mt(it,:,nsim) = m0;
        is_rl = isrl(it);
    end

    if is_rl
        % use RL task ordering

    else
        % use inference task ordering

    end

    
    
end



end