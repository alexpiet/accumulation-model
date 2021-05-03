function [particle] = compute_particles(trial, params, p, a_sign)

% Set up parameters
n = p.n; % number of ensemble particles
numsteps = round(trial.T/p.dt);

% set up variable storage
a = zeros(n, numsteps);
fd = zeros(length(p.avals), numsteps);
pd = zeros(length(p.avals), numsteps);
dd = zeros(length(p.avals), numsteps);

% Run trajectories-------------
% Get adapted clicks
[cl, cr]    = make_adapted_cat_clicks(trial.leftbups, trial.rightbups, params(5), params(6));

%slight timing offset on click times because make_click_inputs35 uses qfind, whereas compute_full_trial rounds click times and then bins. The following code is a little messy because it attempts to counteract qfind. I should really just not use make_click_inputs35, but the error is only off by 1 dt, so it doesn't really matter. 
tvec = 0:p.dt:trial.T;
[difflr, sumlr] = make_click_inputs35(tvec, trial.leftbups, trial.rightbups, cl, cr);
%trial.leftbups = round(trial.leftbups/p.dt)*p.dt;
%trial.rightbups = round(trial.rightbups/p.dt)*p.dt;
%difflr = [difflr(2:end); difflr(end)];
%sumlr = [sumlr(2:end); sumlr(end)];

% Initalize with variance
a(:,1) = sqrt(params(4)).*randn(n,1);

% run foward euler
for i=1:numsteps-1
    a(:,i+1) = a(:,i)...
             + p.dt*params(1).*a(:,i)...
             + difflr(i) + sqrt(sumlr(i)*params(3)).*randn(n,1) ...
             + sqrt(params(2)*p.dt).*randn(n,1);
end

% Compute foward distributions, using hist()
for i=1:numsteps
    counts = hist(a(:,i), p.avals);
    fd(:,i) = counts./n;
end

% compute backward distributions for bins on "a" using a_grid

dDEX = a(:,end)>p.d_dex-0.5 & a(:,end) < p.d_dex+0.5;
nddex = sum(dDEX);
for i=1:numsteps
    counts = hist(a(dDEX,i), p.avals);
    dd(:,i) = counts./nddex;
end

% Compute posterior distributions, using a_sign
if p.a_sign > 0
    DEX = a(:,end)>0;
else
    DEX = a(:,end)<0;
end
ndex = sum(DEX);
for i=1:numsteps
    counts = hist(a(DEX,i), p.avals);
    pd(:,i) = counts./ndex;
end

% pack everything up
particle.fpdf = fd'.*(1/p.da);
particle.ppdf = pd'.*(1/p.da);
particle.dpdf = dd'.*(1/p.da);
particle.a = a;
particle.avals = p.avals;
particle.T = tvec;
particle.DEX = DEX;
particle.dDex = dDEX;
