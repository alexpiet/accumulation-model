function [posterior] = down_sample_posterior(posterior, p)
    ns = size(posterior.ma,2);
    ds_index = p.eval_dt/p.dt;
    if ds_index == 1
        return;
    elseif ds_index < 1
        error('eval dt smaller than dt')
    elseif floor(ds_index) ~= ds_index
        error('eval dt is a non integer multiple of dt')
    end 
    fs = floor(ns/ds_index);    
    good_dex = (0:ds_index:fs*ds_index) + 1;
    if good_dex(end) > ns;
        good_dex(end) = [];
    end   
 
    posterior.ma = posterior.ma(:,good_dex);
    posterior.va = posterior.va(:,good_dex);
    posterior.s  = posterior.s(:,good_dex);
    posterior.T  = posterior.T(:,good_dex);

