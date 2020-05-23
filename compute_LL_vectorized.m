function [NLL_total, mean_a,var_a, NLL, p_chosen, p_chooseR] = ...
   compute_LL_vectorized(buptimes,streamIdx,stim_duration, pokedR, param_set,varargin)

    p=inputParser;
    p.addParameter('nantimes',[]);
    p.addParameter('prior_mean',nan(1,8));
    p.addParameter('prior_var',nan(1,8));
    p.addParameter('use_param', true(1,8));
    p.addParameter('param_default', [0 1 1 0 1 0.1 0 0.01]);
    p.addParameter('adaptation_scales_perclick','var',@(x) any(strcmp(x,{'std','var','none'})));
    parse(p, varargin{:})
    params = p.Results; 


    if isempty(params.nantimes)
        error('you should pass in nantimes')
        nantimes = nan(buptimes);
    else
        nantimes = params.nantimes;
    end
    
    % use default params where necessary
    params_full = params.param_default;
    params_full(logical(params.use_param)) = param_set;    
    
    %% apply adaptation
    if abs(params_full(5) - 1) > eps    
        adapted = adapt_vectorized_clicks(buptimes,nantimes, params_full(5), params_full(6)) .* streamIdx;
    else
        adapted=streamIdx;
        adapted(nantimes)=nan;
    end
    temp = exp(params_full(1)*bsxfun(@minus,stim_duration,buptimes));
    %% compute mean of distribution 
    mean_a = nansum(adapted.*temp);
    %% compute variance of distribution
    % Three components: initial (parameters(4)), accumulation (parameters(2)), and per-click (parameters(3))
    init_var = max(eps,params_full(4));
    a_var = max(eps,params_full(2));
    c_var = max(eps,params_full(3));
    % Initial variance and accumulation variance
    if abs(params_full(1)) < 1e-10
        s2 = init_var*exp(2*params_full(1)*stim_duration) + a_var*stim_duration;
    else
        s2 = init_var*exp(2*params_full(1)*stim_duration) + (a_var./(2*params_full(1)))*(exp(2*params_full(1)*stim_duration)-1);
    end  
    % Add per-click variance
    if strcmp('std',params.adaptation_scales_perclick)
        var_a = s2 + nansum (c_var.*abs(adapted).^2.*temp.^2);
    elseif strcmp('var',params.adaptation_scales_perclick)
        var_a = s2 + nansum (c_var.*abs(adapted).*temp.^2);
    elseif strcmp('none',params.adaptation_scales_perclick)
        var_a = s2 + nansum (c_var.*temp.^2);
    end
    
    %% calulate log likelihood for all trials
    % i see no reason bias can't be allowed to grow beyond -1, +1 
    % bias=min(max(parameters(7),eps-1),1-eps);
    bias = params_full(7);
    
    lapse = min(max(params_full(8),eps),1-eps);
    NLL = zeros(length(pokedR),1);
    erfTerm = erf( (mean_a - bias) ./ sqrt(2*var_a));
    erfTerm(erfTerm==1)=1-eps;
    erfTerm(erfTerm==-1)=eps-1;
    
    
    p_chooseR = ((1-lapse).*(1+erfTerm)+lapse)/2;
    p_chosen = p_chooseR;
    p_chosen(~pokedR) =  (1 - lapse) .* (1 - 0.5 .* (1+erfTerm(~pokedR))) + lapse/2 ;
    NLL = -log(p_chosen);
%     NLL(pokedR) = -log( p_chooseR ) ;
%     NLL(~pokedR) = -log(p_chosen(~pokedR));
    
    prior_cost = 0;
    for pp = 1:length(params.prior_mean)
        if ~isnan(params.prior_mean(pp)) &  ~isnan(params.prior_mean(pp)) 
            prior_cost = prior_cost + (params_full(pp) - params.prior_mean(pp))^2/(2*params.prior_var(pp)^2);
        end
    end
    
    
    if any(~isfinite(NLL))
        error('%g trials had non-finite log likelihoods!',sum(~isfinite(NLL)));
    end
    NLL_total = sum(NLL) + prior_cost;
end