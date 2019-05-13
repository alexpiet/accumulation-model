function [model,p] =accumulation_model(data, params, command, p_in);
% 
% options: fit, evaluate, forward, backward
% Set default parameters
p.dt            = 1e-4; % Timestep for evaluating model
p.error_tolerance = 1e-4;
p.posterior_only = 0;
p.forward_only  = 0;
p.evaluate_only = 0;
p.fit           = 0;
p.compute_dist  = 0;    % 1=compute full pdf
p.eval_tvec     = [];   % Vector of time points to evaluate PDF
%%p.b             = 50; % Will not stay here
p.eval_dt       = 1e-3; % downsample for making PDF
p.just_pdf      = 0;
p.return_backwards = 0;

% use command option to determine what we will do
if nargin > 2
    if strcmp(command, 'posterior'); p.posterior_only = 1;  end;
    if strcmp(command, 'forward');   p.forward_only = 1;    end;
    if strcmp(command, 'evaluate');  p.evaluate_only = 1;   end;
    if strcmp(command, 'fit');       p.fit = 1;             end;
end

% overwrite default parameters with those in p_in
if nargin > 3; 
    names = fieldnames(p_in);
    for i=1:length(names); p = setfield(p, names{i}, getfield(p_in, names{i})); end;
end

% Parse inputs, do checks
if params(3) < 0.5;     error('You might be hitting the low-variance case.'); end;
if length(data) < 1;    error('no data!'); end;
if length(params) ~= 8; error('This model only works with 8 parameters'); end;

if p.fit
    error('I havent implemented the fitting into this master function yet. Check out the function fit_rat_analytical() to fit the model')
end

% for each trial
p_set = p;
for i=1:length(data)
    p = p_set; % resets p.avals, so we use the smallest grid possible for each trial
    if mod(i,100) == 1
        disp(i)
    end
    try 
    % create inputs
    [cl, cr]    = make_adapted_cat_clicks(data(i).leftbups, data(i).rightbups, params(5), params(6));
    data(i).clicks      = [-cl  +cr];
    data(i).times       = [data(i).leftbups data(i).rightbups];
    data(i).times       = round(data(i).times/p.dt)*p.dt;
    data(i).times(data(i).times == 0) = p.dt; % round up clicks that happen in first time window
    data(i).T           = round(data(i).T/p.dt)*p.dt;
    data(i).numsteps    = round(data(i).T/p.dt);
    data(i).dtimes      = round(data(i).times*(1/p.dt));
    
    % Compute forward pass
    [forward] = compute_full_trial(data(i),params,p);
    if p.compute_dist & ~p.posterior_only; 
        forward = compute_pdf(forward, p.avals, p); 
    end;
    
    % determine grid
    p = get_grid(data(i), forward, params, p);
    
    % compute backward pass
    if p.return_backwards
        [model(i).backwards, model(i).posterior]  = compute_backwards_pass(data(i),params,p,forward);
    else
        [~, model(i).posterior]  = compute_backwards_pass(data(i),params,p,forward);
    end
    % Down-Sample Distribution
    model(i).posterior = down_sample_posterior(model(i).posterior,p); 
    
    % Compute Posterior PDF if requested, slow to compute
    if p.compute_dist; 
        model(i).posterior = compute_pdf(model(i).posterior, p.avals,p,'mixture'); 
        model(i).posterior.mean =  mean((model(i).posterior.pdf.*p.da)*repmat(p.avals',1,size(p.avals,2)),2); 
        if p.just_pdf;
            model(i).posterior = rmfield(model(i).posterior, {'ma','va','s','a_grid','grid_s'});
        end
    end;

    % Return Forward model if requested    
    if ~p.posterior_only
        model(i).forward = forward;
    end
    catch
        disp(['failure - ' num2str(i)])
        model(i).posterior = [];
        
    end

end

