function [forward] = compute_full_trial(trial,params,p)
% Computes the forward pass model
% 
% Inputs:
% trial,            a struct in the same format as Bing's model, with fields
%   trial.T         duration of trial
%   trial.leftbups  Times of left clicks
%   trial.rightbups Times of right clicks
%  
% params,           a vector of model parameters
%  
% p,                a struct that contains
% p.dt              time resolution to evaluate the model at
% p.check_final     logical

% Do parameter check
if ~(length(params) == 8)
    error('this function only works for 8 parameter model')
end

% This block got moved to be shared with forward and backward computations
% This function computes both within and across stream adapatation
%[cl, cr]    = make_adapted_cat_clicks(trial.leftbups, trial.rightbups, params(5), params(6));
%clicks      = [-cl  +cr];
%times       = [trial.leftbups trial.rightbups];
%times       = round(times/p.dt)*p.dt;
%trial.T     = round(trial.T/p.dt)*p.dt;
%numsteps    = round(trial.T/p.dt);
%dtimes      = round(times*(1/p.dt));
% after this point, I only need: clicks, dtimes, times, numsteps

% compute mean of distribution
ma = zeros(1,trial.numsteps);
for i=1:length(trial.clicks)
    ind = trial.dtimes(i):length(ma);
    tx  = (trial.times(i):p.dt:trial.T) - trial.times(i);
    this = trial.clicks(i).*exp(params(1).*(tx));
    ma(ind) = ma(ind) + this;
end

% compute variance of distribution
% Three components: initial (params(4)), accumulation (params(2)), and per-click (params(3))
init_var = params(4);
a_var    = params(2);
c_var    = params(3);
va       = zeros(1,trial.numsteps);

% accumulation and initial noise
if abs(params(1)) < 1e-10
    va = init_var.*exp(2*params(1)*(p.dt:p.dt:trial.T)) + a_var.*(p.dt:p.dt:trial.T);
else
    va = init_var.*exp(2*params(1)*(p.dt:p.dt:trial.T)) + (a_var./(2*params(1))).*(exp(2*params(1).*(p.dt:p.dt:trial.T))-1);
end

% Per-click noise
for i=1:length(trial.clicks);
    ind = trial.dtimes(i):length(va);
    tx  = (trial.times(i):p.dt:trial.T) -trial.times(i);
    tc  = trial.clicks(i);
    va(ind) = va(ind) + c_var*abs(tc)*exp(2*params(1)*(tx));
end

% Pack up vector of means, variance, and timepoints. 
forward.ma = ma;
forward.va = va;
forward.T  = p.dt:p.dt:trial.T;

