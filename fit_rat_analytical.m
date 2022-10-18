function [fit, data, savefile] = fit_rat_analytical(ratnames, varargin)
% fit_rat_analytical(ratnames)
% Script for fitting Bing's accumulation model to the dynamic clicks task behavior
% input:
% ratnames      a list of ratnames
%
% searches for mat files with each ratname in the input list in the data directory specified
% by data_dir
%
% if no input is given, loads a synthetic dataset

p = inputParser();
addParameter(p, 'p0',[]);
addParameter(p, 'data_dir','');
addParameter(p, 'results_dir','');
addParameter(p, 'dosave',true);
addParameter(p, 'reload',true);
addParameter(p, 'save_suffix','');
addParameter(p, 'use_param',true(1,8));
addParameter(p, 'param_default',[]);
addParameter(p, 'vectorize',true);
addParameter(p, 'param_names', {'\lambda', '\sigma_a', '\sigma_s', '\sigma_i', ...
                                            '\phi',  '\tau', 'bias', 'lapse'});
addParameter(p,'prior_mean', [nan 0 nan .5    nan nan nan nan]);%[nan 0 nan 0 nan .1059 nan nan]);
addParameter(p,'prior_var',  [nan 5.39 nan 1.87 nan nan nan nan]);%[nan 5.39 nan 1.87 nan .01 nan nan]);
addParameter(p,'use_priors', 1);
addParameter(p, 'lower_bound', [ -40, 0,   0,    0,  0, 0.01, -20, 0]);
addParameter(p, 'upper_bound', [ +40, 200, 5000, 30, 1.2, 5, +20, 1]);
addParameter(p, 'overwrite',0)
addParameter(p, 'TolFun',1e-6)
addParameter(p, 'skipOptimization',false)
parse(p,varargin{:})
par         = p.Results;
data_dir    = par.data_dir;
results_dir = par.results_dir;
save_suffix = par.save_suffix;
skipOptimization = par.skipOptimization;

% Loop over the rats if desired
if iscell(ratnames)
    nrats = length(ratnames);
    for rr = nrats:-1:1
        ratname = ratnames{rr};
        fprintf('Working on %s (%i of %i)...',ratname, nrats-rr+1, nrats)
        fit(rr) = fit_rat_analytical(ratname, varargin{:});
    end
    return
end

% Get the data to fit
if isstruct(ratnames)
    rawdata = ratnames;
    if isfield(rawdata,'ratname')
        ratname = rawdata(1).ratname;
    else
        ratname = ['unknown_' num2str(now)];
    end
    ratnames = ratname; %#ok<NASGU>
    datafile = [];
elseif isnumeric(ratnames) & ratnames < 1
    % negative ratnum loads synthetic dataset for testing purposes
    ratname = ['dataset_' num2str(abs(ratnum))];
    datafile = fullfile(data_dir, ratname);
    data = load(datafile);
elseif ischar(ratnames)
    ratname = ratnames;
    savefile = fullfile(results_dir, ['fit_analytical_' ratname ...
        save_suffix '.mat']);
    if exist(savefile,'file') && ~par.overwrite
        fprintf('loading existing accumulation fit file...')
        load(savefile,'fit');
        return
    end
    datafile = fullfile(data_dir, ratname);
    data = load(datafile);
end
savefile = fullfile(results_dir, ['fit_analytical_' ratname ...
    save_suffix '.mat']);

% set initial parameter values, upper and lower bounds and prior
if isempty(par.p0)
    %p0      = [0.01*randn(1) 50*rand 50*rand    30*rand     rand    .7*rand 0       .5*rand]
    p0      = [0.01*randn(1) rand rand  rand  .5*rand .7*rand 0 .1*rand];
else
    p0 = par.p0;
end

% Look for results and either return or use as starting point for this fit
if exist(savefile,'file')
    if ~par.overwrite
        fprintf('loading existing file...')
        load(savefile,'fit')
        return
    elseif par.reload && isempty(par.p0)
        old_fit = load(savefile);
        p0  = old_fit.fit.final_full;
        fprintf('initializing at last recovered solution')
    end
end

% Determine which parameters to fit and whether to use priors
use_param   = logical(par.use_param);
p0          = p0(use_param);
param_names = par.param_names;
lower_bound = par.lower_bound(use_param);
upper_bound = par.upper_bound(use_param);
if par.use_priors
    prior_mean  = par.prior_mean;
    prior_var   = par.prior_var;
else
    prior_mean  = nan(8,1);
    prior_var   = nan(8,1);
end
extra_args = {'use_param',par.use_param};
if ~isempty(par.param_default)
    extra_args{end+1} = 'param_default';
    extra_args{end+1} = par.param_default;
end

% Setup the data and the function to minimize
rawdata     = data.rawdata;
pokedR      = [rawdata.pokedR]';
stim_dur    = [rawdata.T];
if par.vectorize
    [buptimes,nantimes,streamIdx] = vectorize_clicks({rawdata.leftbups}, {rawdata.rightbups});
    nllfun = @(params) compute_LL_vectorized(buptimes,streamIdx,...
        stim_dur, pokedR, params,...
        'prior_var',prior_var,'prior_mean',prior_mean,...
        'nantimes', nantimes, extra_args{:});
else
    opts.compute_full = false;
    opts.dt = .001;
    fill_params = @(params) fill_defaults(params,par.use_param,...
        par.param_default);
    nllfun = @(params)  compute_LL(rawdata, fill_params(params), opts);
end
bf_full = p0;

% Find the best fit parameters 
tfit = nan;
if ~skipOptimization
    disp(['Starting fit for rat ' ratname])
    tic
    [xbf, f, exitflag, output, ~, grad, hessian] = ...
        fmincon(nllfun, p0, ...
        [], [], [], [], lower_bound,  upper_bound, [], ...
        optimset('Display', 'iter-detailed','TolFun',par.TolFun,...
        'MaxFunEvals',1e4));
    bf_full(use_param) = xbf;
    tfit = toc;
else
    xbf = bf_full(use_param);
    f = nan;
    exitflag = '';
    grad = [];
    hessian = [];
end


% compute values at best fit solution
if par.vectorize
    [~, ma, va, ~, ~, pr] = compute_LL_vectorized(buptimes,streamIdx,...
        stim_dur, pokedR, xbf,...
        'prior_var',prior_var,'prior_mean',prior_mean,...
        'nantimes', nantimes, extra_args{:});
else
    [~, ma, va, pr] = compute_LL(rawdata, fill_params(xbf_full), opts);
end

% Package outputs
fit.p0          = p0;
fit.lower_bound = lower_bound;
fit.upper_bound = upper_bound;
fit.rat         = ratname;
fit.time        = tfit;
fit.final       = xbf;
fit.f           = f;
fit.final_full  = bf_full;
fit.exitflag    = exitflag;
fit.grad        = grad;
fit.hessian     = hessian;
fit.nt          = length(pokedR);
fit.fpt         = exp(-fit.f/fit.nt);
fit.pr          = pr;
fit.a_mu        = ma;
fit.a_var       = va;
fit.datafile    = datafile;
fit.datenum     = datetime;
fit.use_param   = par.use_param;
fit.param_default = par.param_default;
fit.param_names   = param_names;
fit.prior_mean      = prior_mean;
fit.prior_var      = prior_var;
fit.nllfun      = nllfun;

% Save results
if par.dosave && ~skipOptimization
    save(savefile,'fit','par','data');
end


function param_set = fill_defaults(params,use_param,defaults)
param_set  = defaults;
param_set(logical(use_param)) = params;
