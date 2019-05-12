function [mean_a, var_a] = compute_backwards_trajectory(trial,params,p,T,init_a)
% computes single delta-mode backwards/posterior solution

nd = length(init_a);

% Compute mean of distribution
% (1) displacement from initial conditions
mean_a = init_a'.*exp(params(1).*(T - trial.T));

% (2) displacement from each click
%for i=length(trial.clicks):-1:1
%    mean_a(:,trial.sorted_dtimes(i)) = mean_a(:,trial.sorted_dtimes(i)) - trial.sorted_clicks(i);
%    mean_a(:,1:trial.sorted_dtimes(i)) = mean_a(:,trial.sorted_dtimes(i)).*exp(params(1).*(T(1:trial.sorted_dtimes(i)) - trial.sorted_times(i)));
%end
for i=length(trial.clicks):-1:2
    mean_a(:,trial.sorted_dtimes(i)) = mean_a(:,trial.sorted_dtimes(i)) - trial.sorted_clicks(i);
    mean_a(:,trial.sorted_dtimes(i-1):trial.sorted_dtimes(i)) = mean_a(:,trial.sorted_dtimes(i)).*exp(params(1).*(T(trial.sorted_dtimes(i-1):trial.sorted_dtimes(i)) - trial.sorted_times(i)));
end
mean_a(:,trial.sorted_dtimes(1)) = mean_a(:,trial.sorted_dtimes(1)) - trial.sorted_clicks(1);
mean_a(:,1:trial.sorted_dtimes(1)) = mean_a(:,trial.sorted_dtimes(1)).*exp(params(1).*(T(1:trial.sorted_dtimes(1)) - trial.sorted_times(1)));

% Compute variance of distribution
% Three components: initial (params(4)), accumulation (params(2)), and per-click (params(3))
init_var = params(4);
a_var    = params(2);
c_var    = params(3);
var_a    = zeros(size(mean_a));

% (1) Initial variance and (2) Accumulation variance
if abs(params(1)) < 1e-10
    var_a = repmat(a_var,nd,1).*(trial.T - T);
else
    var_a = repmat(a_var./(2*params(1)),nd,1).*(exp(2*params(1).*(trial.T-T))-1);
end

% (3) Add per-click variance
for i=length(trial.clicks):-1:1
    var_a(:,1:trial.sorted_dtimes(i)) = var_a(:,1:trial.sorted_dtimes(i)) + c_var*abs(trial.sorted_clicks(i))*exp(2*params(1).*(T(1:trial.sorted_dtimes(i)) - trial.sorted_times(i)));
end
 
% add initial variance
var_a(:,1) = var_a(:,1) + init_var;

% force variance to be positive. 0 variance breaks shitty gmdistribution.
var_a(var_a <= 0) = eps;




