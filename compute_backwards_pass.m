function [back, posterior] = compute_backwards_pass(trial,params,p, forward)

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
%times = times - p.dt;
%trial.T     = round(trial.T/p.dt)*p.dt;
%numsteps    = round(trial.T/p.dt);
%dtimes      = round(times*(1/p.dt));
%trial.times = times;
%trial.dtimes = dtimes;
%trial.clicks = clicks;

%%%%% Shifting one earlier for alignment
% We have to do this because we don't want the increase in variance for the forward and backward models 
% to kick in on the same timestep, because then we'd get double variance on that timestep. So we pick
% one to offset, and I picked the backwards model. 
trial.times  = trial.times - p.dt;
trial.dtimes = trial.dtimes - 1;

% remove trials that get pushed to 0
rdex = trial.dtimes == 0;
trial.times(rdex) = [];
trial.dtimes(rdex) = [];
trial.clicks(rdex) = [];

% sorted times are useful for iterative solution construction
[trial.sorted_times, trial.sort_index] =sort(trial.times);
trial.sorted_dtimes = trial.dtimes(trial.sort_index);
trial.sorted_clicks = trial.clicks(trial.sort_index);

% set up some discritization of the decision axis. 
if trial.pokedR > 0
    init_as = p.a_grid(p.a_grid >= params(7));
    grid_s = p.a_grid_s(p.a_grid >=params(7));
else
    init_as = p.a_grid(p.a_grid <= params(7));
    grid_s = p.a_grid_s(p.a_grid<= params(7));
end

% Compute delta-solution at each grid point
%for i=1:length(init_as)
%    [back.ma(i,:), back.va(i,:)] = compute_backwards_trajectory_single(trial, params, p,forward.T, init_as(i));
% iterative solution here for reference while things get vectorized
%    posterior.ma(i,:) = (forward.ma.*back.va(i,:) + back.ma(i,:).*forward.va)./(forward.va + back.va(i,:));
%    posterior.va(i,:) = (forward.va.*back.va(i,:))./(forward.va + back.va(i,:));
%    posterior.s(i,:) = (1./sqrt(2*pi*(forward.va + back.va(i,:)))).*exp(-(forward.ma - back.ma(i,:)  ).^2./(2.*(forward.va+back.va(i,:))));
%    posterior.s(i,:) = posterior.s(i,:).*grid_s(i); 
%end

% Compute delta-solution at each grid point

[back.ma, back.va] = compute_backwards_trajectory(trial, params, p,forward.T, init_as);

% Compute posterior distribution for each delta-solution function
% Need to compute mean, variance, and scale-factor for each product of forward*backward-delta-solution
% scale based on grid density and alignment
% Equation for scale factor for reference:
% scale = 1/(sqrt(2*pi*(s2f+s2g)))*exp(-(muf-mug)^2/(2*(s2f+s2g))
posterior.ma = (forward.ma.*back.va + back.ma.*forward.va)./(forward.va + back.va); 
posterior.va = (forward.va.*back.va)./(forward.va + back.va);
posterior.s  = (1./sqrt(2*pi*(forward.va + back.va))).*exp(-(forward.ma - back.ma  ).^2./(2.*(forward.va+back.va)));
posterior.s  = posterior.s.*repmat(grid_s', 1, size(posterior.s,2));

% Put the 0 mode at half strength. 
zero_dex = find(init_as == params(7));
if ~isempty(zero_dex) 
    if length(zero_dex) > 1
        % Throw an error if we have multiple 0 modes. 
        error('redundant zero mode!')
    end
    posterior.s(zero_dex,:) = 0.5.*posterior.s(zero_dex,:);
end

% force scaling weights to be > 0, shitty gmdistribution can't handle ==0
posterior.s(posterior.s <= 0) = eps;

% Pack up discretization for use later
back.T      = forward.T;
back.a_grid = init_as;
back.grid_s = grid_s;
posterior.T = forward.T;
posterior.a_grid = init_as;
posterior.grid_s = grid_s;


