function [] = fit_rat_analytical(ratnum)
% Script for fitting Bing's accumulation model to the dynamic clicks task behavior

rats = {'high_gamma/H065','high_gamma/H067'};

% Try to fit model
try
%    parpool
    % Load data
    if ratnum < 1
        % negative ratnum loads synthetic dataset for testing purposes
        ratname = ['dataset_' num2str(abs(ratnum))];
        load(['../DATA/' ratname]);
    else
        ratname = rats{ratnum};
        load(['../DATA/' ratname]);% '_' num2str(dex)])
    end
    
    % set all other parameters
    params      = [0.01*randn(1) 50*rand 50*rand    30*rand     rand    .7*rand 0       .5*rand];
    %params =     [ lambda,     sa,     ss,         si,         phi,    tau,    bias,   lapse];
    lower_bound = [ -40,       0,      0,           0,          0,      0,      -.1,     0];
    upper_bound = [ +40,       200,    5000,        30,         1.2,    0.7,    +.1,     1];
    p.prior     = [0,           5.39,   0,          1.87,       0,      0,      0,      0];
    
    p.compute_full = 0;
    p.dt = 1e-4;

    % fit
    disp(['Starting fit for rat ' ratname])
    tic
    [x_fmincon, f, exitflag, output, ~, grad, hessian] = fmincon(@(parameters) compute_LL(data, parameters, p), params, [], [], [], [], lower_bound,  upper_bound, [], optimset('Display', 'iter-detailed'));
    t=toc;

    fit.rat = ratname;
    fit.time = t;
    fit.final = x_fmincon;
    fit.f = f;
    fit.exitflag = exitflag;
    fit.grad = grad;
    fit.hessian = hessian;

    % Save results
    savefile = ['../OUTPUT/fit_analytical_' ratname ];
    save(savefile);
catch me
    disp(me)

    fprintf('Rat failed: %s\n',ratname)
end


