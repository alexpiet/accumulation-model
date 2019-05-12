% Settings

ratname = 'H033';
param_num = 4;
param_vals = [0:10:200];

%%%
load(['../FITS/fit_analytical_' ratname])
params = fit.final;
NLLs = zeros(1,length(param_vals));

for i=1:length(param_vals)
    params(param_num) = param_vals(i);
    NLLs(i) = compute_LL(data,params,p);
end

dNLL = NLLs - fit.f;


