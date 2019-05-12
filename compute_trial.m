function [mean_a, var_a] = compute_trial(trial,params,p)

% Need these things:
%trial.pokedR
%trial.leftbups
%trial.right
%trial.T

%params = lambda, sa, ss, si, phi, tau, bias, lapse

%p.dt
%p.compute_full = 1 (entire time trajectory), 2 (just end value)


% Compute list of effective clicks at each times by computing adaptation

% This function only computes within stream adaptation
%[cl, cr] = make_adapted_clicks(trial.leftbups, trial.rightbups, params(5), params(6), 1);

% This function computes both within and across stream adapatation
if length(params) == 8
    [cl, cr] = make_adapted_cat_clicks(trial.leftbups, trial.rightbups, params(5), params(6));
elseif length(params) == 7
    [cl, cr] = make_adapted_cat_clicks(trial.leftbups, trial.rightbups, params(4), params(5));
elseif length(params) == 6
    [cl, cr] = make_adapted_cat_clicks(trial.leftbups, trial.rightbups, params(3), params(4));
else
    error('weird number of parameters')
end

clicks = [-cl  +cr];
times  = [trial.leftbups trial.rightbups];

% compute mean of distribution
if p.compute_full
    error('not implemented')
    numsteps = round(trial.T/p.dt);
    mean_a = zeros(1,numsteps);
    ts = 0;...; % time of each click in timestep indexing
    for i=1:length(clicks)
        mean_a(ts(i):end) = mean_a(ts(i):end) + clicks(i).*exp(params(1).*((times(i):p.dt:trial.T) - times(i)));
    end
else
    mean_a = 0;
    for i=1:length(clicks)
        mean_a = mean_a + clicks(i)*exp(params(1)*(trial.T - times(i)));
    end
end


% compute variance of distribution
% Three components: initial (params(4)), accumulation (params(2)), and per-click (params(3))
if p.compute_full
    error('not implemented')
else
    if length(params) == 8
        init_var = params(4);
    else
        init_var = 0;
    end
    if length(params) == 6
        a_var = 0;
        c_var = params(2);
    else
        a_var = params(2);
        c_var = params(3);
    end


    % Initial variance and accumulation variance
    if abs(params(1)) < 1e-10
        s2 = init_var*exp(2*params(1)*trial.T) + a_var*trial.T;
    else
        s2 = init_var*exp(2*params(1)*trial.T) + (a_var./(2*params(1)))*(exp(2*params(1)*trial.T)-1);
    end
    
    % Add per-click variance
    for i=1:length(clicks)
        s2 = s2 + c_var*abs(clicks(i))*exp(2*params(1)*(trial.T - times(i)));
    end
end
var_a = s2;

