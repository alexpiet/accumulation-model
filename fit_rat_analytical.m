function [] = fit_rat_analytical(ratnames)
% fit_rat_analytical(ratnames)
% Script for fitting Bing's accumulation model to the dynamic clicks task behavior
% input:
% ratnames      a list of ratnames 
%
% searches for mat files with each ratname in the input list in the data directory specified
% by pbups_dyn_path
%
% if no input is given, loads a synthetic dataset


pbups_dyn_path()

if iscell(ratnames) 
    for rr = 1:length(ratnames)
        ratname = ratnames{rr};
        fit_rat_analyical(ratname);
    end
else
    
    
    
    % Try to fit model
    try
        %    parpool
        % Load data
        if isnumeric(ratnames) & ratnames < 1
            % negative ratnum loads synthetic dataset for testing purposes
            ratname = ['dataset_' num2str(abs(ratnum))];
            load(fullfile(data_dir, ratname));
        else
            ratname = ratnames;
            if exist(fullfile(example_data_dir, [ratname '.mat']),'file')
                load(fullfile(example_data_dir, ratname));
            elseif exist(fullfile(data_dir, [ratname '.mat']),'file')
                load(fullfile(data_dir, ratname));
            end
        end
        
        % set initial parameter values, upper and lower bounds and prior
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
        [x_fmincon, f, exitflag, output, ~, grad, hessian] = ...
            fmincon(@(parameters) compute_LL(data, parameters, p), params, ...
            [], [], [], [], lower_bound,  upper_bound, [], optimset('Display', 'iter-detailed'));
        t=toc;
        
        fit.rat = ratname;
        fit.time = t;
        fit.final = x_fmincon;
        fit.f = f;
        fit.exitflag = exitflag;
        fit.grad = grad;
        fit.hessian = hessian;
        
        % Save results
        savefile = fullfile(results_dir, ['fit_analytical_' ratname]);
        save(savefile);
    catch me
        disp(me)
        
        fprintf('Rat failed: %s\n',ratname)
    end
end

