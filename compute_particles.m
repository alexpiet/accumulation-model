function [particle] = compute_particles(trial, params, p)

% Set up parameters
n = p.n; % number of ensemble particles
numsteps = round(trial.T/p.dt);

% set up variable storage
a = zeros(n, numsteps);
b = zeros(n, numsteps);
fd = zeros(length(p.avals), numsteps);
pd = zeros(length(p.avals), numsteps);
bd = zeros(length(p.avals), numsteps);
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

% compute backward-delta distribution for bins on "a" using a_grid
%% NEED TO UPDATE HERE. With non-uniform grid, this might be buggy
dDEX = a(:,end)>p.d_dex-p.da_grid/2 & a(:,end) < p.d_dex+p.da_grid/2;
nddex = sum(dDEX);
for i=1:numsteps
    counts = hist(a(dDEX,i), p.avals);
    dd(:,i) = counts./nddex;
end

% Compute posterior distributions, using a_sign
if trial.pokedR > 0
    DEX = a(:,end)>params(7);
else
    DEX = a(:,end)<params(7);
end
ndex = sum(DEX);
for i=1:numsteps
    counts = hist(a(DEX,i), p.avals);
    pd(:,i) = counts./ndex;
end

% Have to shift times off by 1 dt to deal with forward/backwards issue
% but, because times are off by 1 dt for particles anyways, it doesnt matter
%difflr2 = [difflr(2:end); 0];
%sumlr2 = [sumlr(2:end); 0];

% Compute backwards particle trajectories
% Initialize randomly but uniformly across grid bin
%% NEED TO UPDATE HERE. With non uniform grid, this might be buggy
b(:,end) = (p.d_dex - p.da_grid/2) + p.da_grid.*rand(n,1);
for i=numsteps:-1:2
    b(:,i-1) =  b(:,i) ...
                + p.dt*(-params(1)).*b(:,i)...
                - difflr(i-1) + sqrt(sumlr(i-1)*params(3)).*randn(n,1) ...
                + sqrt(params(2)*p.dt).*randn(n,1);
end

% Need to add initial noise on last time step
b(:,1) = b(:,1) + sqrt(params(4)).*randn(n,1);

% compute backwards-delta-particle distribution
for i=1:numsteps
    counts = hist(b(:,i),p.avals);
    bd(:,i) = counts./n;
end


% pack everything up
particle.fpdf = fd'.*(1/p.da);
particle.ppdf = pd'.*(1/p.da);
particle.dpdf = dd'.*(1/p.da);
particle.bpdf = bd'.*(1/p.da);
particle.a = a;
particle.b = b;
particle.avals = p.avals;
particle.T = tvec;
particle.DEX = DEX;
particle.dDex = dDEX;
