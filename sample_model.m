function [data] = sample_model(data,fit,p);

bias =  fit.final(7);
lapse = fit.final(8);

for i=1:length(data)
    [ma,va] = compute_trial(data(i), fit.final,p);
    
    % compute pr, pl with bias
    pr = 0.5*(1+erf( -(bias-ma)/sqrt(2*va)));
%%    pl = 1-pr;

    % compute pr, pl with lapse
    PR = (1-lapse)*pr + lapse*0.5;
%%    PL = (1-lapse)*pl + lapse*0.5;

    data(i).pokedR = rand < PR;
end


