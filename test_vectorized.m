% load the old fit for this rat
rat = 'H084';
f = load('fit_analysis_analyticalH084.mat')
saved_bf = f.fit.final
%% refit with vectorized analytical model
r = fit_rat_analytical(rat)
%%
r = load('../../results/fit_analytical_H084.mat');
vector_bf = r.fit.final;
%%
d = load('H084.mat') ;
data = d.data;
pokedR      = [data.pokedR]';
stim_dur    = [data.T]';
%%
p.dt = 1e-4;
p.compute_full = 0;
prior_mean  = nan(1,8);
prior_var   = nan(1,8);
prior_mean([2 4]) = [0 0];
prior_var([2 4]) = [ 0.01 0.01];
p.prior = nan(1,8);
p.prior([2 4]) = prior_var([2 4]);

test_params = saved_bf;

[nll cost] = compute_LL(data,test_params,p);


[vectorized_nll,~,~,~,vectorized_cost] = compute_LL_vectorized(buptimes,streamIdx,...
    stim_dur', pokedR, test_params, 'nantimes', nantimes,...
    'prior_mean',prior_mean, 'prior_var',prior_var);

nll 
vectorized_nll
cost
vectorized_cost


%%
%params =     [ lambda,     sa,     ss,         si,         phi,    tau,    bias,   lapse];


p.dt = 1e-4;
p.prior = [0,          5.39,   0,          1.87,       0,      0,      0,      0];
p.prior = [0,          5.39,   0,          1.87,       0,      0,      0,      0];
[buptimes,nantimes,streamIdx] = vectorize_clicks({data.leftbups}, {data.rightbups});



buptimes = buptimes;
streamIdx = streamIdx;
nantimes = nantimes;

%%
params       = f.fit.final;
%params = [0 0 0 0 .5 0.01 0 0];

ii = randi(length(data));

[mean_a, var_a] = compute_trial(data(ii), params, p);
[NLL_total,vectorized_mean_a,vectorized_var_a,NLL] = compute_LL_vectorized(...
    buptimes(:,ii),streamIdx(:,ii),...
    stim_dur(ii), pokedR(ii), params, 'nantimes', nantimes(:,ii));

[mean_a, vectorized_mean_a]
[var_a, vectorized_var_a]
%%
try
p = rmfield(p, 'prior')
catch
end

nll = compute_LL(data,params,p);

[vectorized_nll] = compute_LL_vectorized(buptimes,streamIdx,...
    stim_dur', pokedR, params, 'nantimes', nantimes);

[nll vectorized_nll nll-vectorized_nll]

%%


%% compare adaptations
ii = 1;
phi = .5;
tau_phi = .01;
trial = data(ii);
[cl, cr] = make_adapted_cat_clicks(trial.leftbups, trial.rightbups, phi, tau_phi);

adapted = adapt_vectorized_clicks(buptimes(ii,:)',nantimes(ii,:)', phi, tau_phi);
cl_ = abs(adapted(streamIdx(ii,:)==-1));
cr_ = adapted(streamIdx(ii,:)==1);
[sum(cl_) - sum(cr_), sum(cl) - sum(cr)]
%%
sum(cl_) == sum(cl)

%%

%%


adapted = adapt_vectorized_clicks(buptimes,nantimes, params(5), params(6));
whos adapted

%%

%%
[mean_a, vectorized_mean_a]

[var_a, vectorized_var_a]

