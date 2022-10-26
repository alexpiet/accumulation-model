function [NLL_total, ma, va,  p_chooseR, NLL] = compute_LL(data, params,p)


NLL   = zeros(1,length(data));
if length(params) == 8
    bias  = params(7);
    lapse = params(8);
elseif length(params) == 7
    bias = params(6);
    lapse = params(7); 
elseif length(params) == 6
    bias = params(5);
    lapse = params(6);
end

% iterate over trials
p_chooseR = nan(length(data),1);
ma = nan(size(p_chooseR));
va = nan(size(p_chooseR));
for tt=1:length(data)
    [ma(tt),va(tt)] = compute_trial(data(tt), params,p);
    
    % compute pr, pl with bias
    pr = 0.5*(1+erf( -(bias-ma(tt))/sqrt(2*va(tt))));
    
    p_chooseR(tt) = pr;
    
    if pr == 1
        pr = pr - eps;
    elseif pr == 0
        pr = pr + eps;
    end
    
    pl = 1-pr;

    % compute pr, pl with lapse
    PR = (1-lapse)*pr + lapse*0.5;
    PL = (1-lapse)*pl + lapse*0.5;

    % compute NLL for this trial
    if data(tt).pokedR
        nll = -log(PR);
    else
        nll = -log(PL);
    end
    
    % add to total over all trials
    NLL(tt) = nll;
end

NLL_total = sum(NLL);

if isfield(p, 'prior')
    cost = params(2)^2/(2*p.prior(2)^2) + params(4)^2/(2*p.prior(4)^2);
%cost = params(2)^2/(2*p.prior(2)^2) + params(4)^2/(2*p.prior(4)^2) + (params(5) - p.mean_prior(5))^2/(2*p.prior(5)^2) + (params(6)-p.mean_prior(6))^2/(2*p.prior(6)^2) + (params(3)-p.mean_prior(3))^2/(2*p.prior(3)^2)  ;

    NLL_total = NLL_total + cost;
end



