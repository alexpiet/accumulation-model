function [fit] = fit_rat_analytical(ratnames, dosave)
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
        elseif isstr(ratnames)
            ratname = ratnames;
            if exist(fullfile(example_data_dir, [ratname '.mat']),'file')
                load(fullfile(example_data_dir, ratname));
            elseif exist(fullfile(data_dir, [ratname '.mat']),'file')
                load(fullfile(data_dir, ratname));
            end
        elseif isstruct(ratnames)
            data = ratnames;
            ratname = 'unknown';
        end
        
        % set initial parameter values, upper and lower bounds and prior
        p0      = [0.01*randn(1) 50*rand 50*rand    30*rand     rand    .7*rand 0       .5*rand];
        %params =     [ lambda,     sa,     ss,         si,         phi,    tau,    bias,   lapse];
        lower_bound = [ -40,       0,      0,           0,          0,      0,      -20,     0];
        upper_bound = [ +40,       200,    5000,        30,         1.2,    0.7,    +20,     1];
        prior_mean  = nan(1,8);
        prior_var   = nan(1,8);
        prior_mean([2 4]) = [0 0];
        prior_var([2 4]) = [5.39 1.87];
%         p.compute_full = 0;
%         p.dt = 1e-4;
%         
        % fit
        disp(['Starting fit for rat ' ratname])
        tic
        
        [buptimes,nantimes,streamIdx] = vectorize_clicks({data.leftbups}, {data.rightbups});
        pokedR      = [data.pokedR]';
        stim_dur    = [data.T];
        nllfun = @(params) compute_LL_vectorized(buptimes,streamIdx,stim_dur, pokedR, params,...
            'prior_var',prior_var,'prior_mean',prior_mean, 'nantimes', nantimes);
        
        
        [x_fmincon, f, exitflag, output, ~, grad, hessian] = ...
            fmincon(nllfun, p0, ...
            [], [], [], [], lower_bound,  upper_bound, [], ...
            optimset('Display', 'iter-detailed','TolFun',1e-12));
        t=toc;
        
        fit.p0          = p0;
        fit.lower_bound = lower_bound;
        fit.upper_bound = upper_bound;
        fit.rat         = ratname;
        fit.time        = t;
        fit.final       = x_fmincon;
        fit.f           = f;
        fit.exitflag    = exitflag;
        fit.grad        = grad;
        fit.hessian     = hessian;
        fit.nt          = length(pokedR);
        
        % Save results
        savefile = fullfile(results_dir, ['fit_analytical_' ratname]);
        if dosave
            save(savefile);
        end
    catch me
        disp(me)
        
        fprintf('Rat failed: %s\n',ratname)
    end
end

