"""
    compute_LL(data, params; prior_mean=[], prior_var=[])

    Compute the negative log likelihood of the model given the data and parameters

"""
function compute_LL(data, params; prior_mean=[], prior_var=[])
# Set up variables
NLL = 0;
if length(params) == 8
    bias  = params[7];
    lapse = params[8];
elseif length(params) == 7
    bias = params[7];
    lapse = 0;
else
    bias  = params[5];
    lapse = params[6];
end

# iterate over trials
nt = length(data["pokedR"])
for tt=1:nt
    ma, va = compute_trial(data,tt,params);
    # compute pr, pl with bias
    this_pr = 0.5*(1+erf( -(bias-ma)/sqrt(2*va)));
    this_pl = 1-this_pr;

    # compute pr, pl with lapse
    PR = (1-lapse)*this_pr + lapse*0.5;
    PL = (1-lapse)*this_pl + lapse*0.5;
    
    # compute NLL for this trial
    if data["pokedR"][tt] 
        nll = -log(PR);
    else
        nll = -log(PL);
    end
    NLL += nll
end

# add prior cost
if !isempty(prior_mean)
    NLL += compute_prior_cost(params, prior_mean, prior_var)
end

# Return NLL
return NLL
end

function compute_prior_cost(params, prior_mean, prior_var)
    prior_cost = 0
    for pp = 1:length(prior_mean)
        if !isnan(prior_mean[pp])
            prior_cost += (params[pp]-prior_mean[pp])^2 / (2*prior_var[pp]^2)
        end
    end
    return prior_cost
end


"""
    compute_LL(ratname; res_dir="./", data_dir="./", overwrite=true)

    Compute the negative log likelihood of the model given the data and parameters

"""
function compute_LL(ratname; res_dir="./", data_dir="./", overwrite=true)


    # load rat data
    println(res_dir)
    println(ratname)
    fitfn = joinpath(res_dir,"fit_analytical_"*ratname*".mat");
    # load fit if it exists
    if isfile(fitfn)
        fit_data = matread(fitfn);
        fit      = fit_data["fit"] 
    else
        error("couldn't find fit file "*fitfn)
    end
    println("fit loaded")
    fit     = fit_data["fit"];
    println("found field fit")
    data    = fit_data["data"]
    if haskey(data,"rawdata")
        data = data["rawdata"]
    end
    println("found field data")
    params  = fit["final"];
    println("found field final")
    prior_mean = fit["prior_mean"]
    prior_var = fit["prior_var"]
    println("found prior mean and variance")
    nt = length(data["pokedR"])
    println(string(nt)*" trials in this dataset")
    # evaluate model LL just to make sure its correct
    NLL, res = compute_LL_res(data, params, prior_mean = prior_mean, prior_var = prior_var);
    println(NLL)
    bad_NLL = abs(NLL -  fit["f"]) > 1
    if bad_NLL
        println("Oh no! Julia version is not within tolerance of MATLAB values")
        temp = NLL - fit["f"];
        println(temp)
    else
        println("Good. Julia NLL is within tolerance of MATLAB values")
    end
    return res
end


"""
    compute_LL_res(data, params; prior_mean=[], prior_var=[])

    Compute the negative log likelihood of the model given the data and parameters
"""
function compute_LL_res(data, params; prior_mean=[], prior_var=[])
# Set up variables
NLL = 0;
if length(params) == 8
    bias  = params[7];
    lapse = params[8];
elseif length(params) == 7
    bias = params[7];
    lapse = 0;
else
    bias  = params[5];
    lapse = params[6];
end

# iterate over trials
nt = length(data["pokedR"])
a_mu = Array{Float64, 1}(undef, nt)
a_var = Array{Float64, 1}(undef, nt)
nll = Array{Float64, 1}(undef, nt)
pr = Array{Float64, 1}(undef, nt)

for tt=1:nt
    ma, va = compute_trial(data,tt,params);
    # compute pr, pl with bias
    this_pr = 0.5*(1+erf( -(bias-ma)/sqrt(2*va)));
    this_pl = 1-this_pr;

    # compute pr, pl with lapse
    PR = (1-lapse)*this_pr + lapse*0.5;
    PL = (1-lapse)*this_pl + lapse*0.5;

    PR = min(1-eps(), max(eps(), PR))
    PL = min(1-eps(), max(eps(), PL))
    
    # compute NLL for this trial
    if data["pokedR"][tt] 
        nll[tt] = -log(PR);
    else
        nll[tt] = -log(PL);
    end
    a_mu[tt] = ma
    a_var[tt] = va
    pr[tt] = PR
    
    NLL += nll[tt]
end

# add priors
prior_cost = 0
for pp = 1:length(prior_mean)
    if !isnan(prior_mean[pp])
        prior_cost += (params[pp]-prior_mean[pp])^2 / (2*prior_var[pp]^2)
    end
end
NLL += prior_cost

res =  Dict("a_mu"=>a_mu, "a_var"=>a_var, "pr"=>pr, "nll"=>nll, "NLL"=>NLL)
return NLL, res
end





"""
    compute_trial(data, i, params)

    Compute the mean and variance of the distribution of the accumulation process at the time of choice for a given trial
"""
function compute_trial(data, i, params);
    
    # run clicks through the adaptation process  
    if length(params) == 8
        cl, cr = make_adapted_cat_clicks(data["leftbups"][i], data["rightbups"][i], params[5],params[6]);
    elseif length(params) == 7
#        cl, cr = make_adapted_cat_clicks(data["leftbups"][i], data["rightbups"][i], params[4],params[5]);
        cl, cr = make_adapted_cat_clicks(data["leftbups"][i], data["rightbups"][i], params[5],params[6]);
    elseif length(params) == 6;
        cl, cr = make_adapted_cat_clicks(data["leftbups"][i], data["rightbups"][i], params[3],params[4]);
    end

    if !isempty(cl) & !isempty(cr)    
        clicks = [-cl cr];
        times = [data["leftbups"][i] data["rightbups"][i]];
    elseif isempty(cl) & !isempty(cr)
        clicks = cr;
        times = data["rightbups"][i];
    elseif !isempty(cl) & isempty(cr)
        clicks = -cl;
        times = data["leftbups"][i];
    else
        clicks = [];
        times = [];
    end

    # compute mean of distribution
    mean_a = 0;
    for j=1:length(clicks)
        try
            mean_a += clicks[j]*exp(params[1]*(data["T"][i]-times[j]));
        catch
            println("Error in compute_trial")
            println("clicks:"*string(clicks))
            println("times:"*string(times))
            println("T:"*string(data["T"][i]))
            println("lambda:"*string(params[1]))
            println("time j:"*string(times[j]))
        end
    end
    # compute variance of distribution
    # three sources: initial (params[4]), accumulation (params[2]), and per-click (params[3])
    
    if length(params) == 8
        a_var    = params[2];
        c_var    = params[3];
        init_var = params[4];
    elseif length(params) == 7
#        a_var    = params[2];
#        c_var    = params[3];
#        init_var = 0;
        a_var    = params[2];
        c_var    = params[3];
        init_var = params[4];
    elseif length(params) == 6
        a_var    = 0;;
        c_var    = params[2];
        init_var = 0;
    end
    
    # Initial and accumulation variance
    if abs(params[1]) < 1e-10
        s2 = init_var*exp(2*params[1]*data["T"][i]) + a_var*data["T"][i];
    else
        s2 = init_var*exp(2*params[1]*data["T"][i]) + (a_var/(2*params[1]))*(exp(2*params[1]*data["T"][i])-1);
    end
    
    # add per-click variance
    for j=1:length(clicks)
        s2 += c_var*abs(clicks[j])*exp(2*params[1]*(data["T"][i] - times[j]));
    end

    var_a = s2;

    # return mean and variance of distribution
    return mean_a, var_a
end



"""
    make_adapted_clicks(leftbups, rightbups, phi, tau_phi, psi, tau_psi)

    Adaptation function with separate within and across stream adaptation

"""
function make_adapted_clicks(leftbups, rightbups, phi, tau_phi, psi, tau_psi)
    Lsame = ones(typeof(phi),size(leftbups));
    Rsame = ones(typeof(phi),size(rightbups));
    """
    # magnitude of stereo clicks set to zero
    if ~isempty(leftbups) && ~isempty(rightbups) && abs(leftbups[1]-rightbups[1]) < eps()
        Lsame[1] = 0;
        Rsame[1] = 0;
    end;
    """
    # if there's appreciable same-side adaptation
    if abs(phi - 1) > eps() 
        # inter-click-intervals
    #    ici_L = diff(leftbups')';
    #    ici_R = diff(rightbups')';
        if length(leftbups) <= 1
            ici_l = [];
        else
            ici_L = (leftbups[2:end]  - leftbups[1:end-1])';
        end
    
        if length(rightbups) <= 1
            ici_R = [];
        else
            ici_R = (rightbups[2:end]  - rightbups[1:end-1])';
        end
        
        for i = 2:length(leftbups),
            last_L = tau_phi*log(1-Lsame[i-1]*phi);
            Lsame[i] = 1 - exp((-ici_L[i-1] + last_L)/tau_phi);
        end;
        
        for i = 2:length(rightbups),
            last_R = tau_phi*log(1-Rsame[i-1]*phi);
            Rsame[i] = 1 - exp((-ici_R[i-1] + last_R)/tau_phi);
        end;
        
        Lsame = real(Lsame);
        Rsame = real(Rsame);
    end;
    
    Lother = ones(size(leftbups));
    Rother = ones(size(rightbups));
    
    # if there's appreciable across-side adaptation
    if abs(psi - 1) > eps() 
    #    lefts  = [leftbups(:)  -ones(numel(leftbups),1)];
    #    rights = [rightbups(:) +ones(numel(rightbups),1)];
    #    allbups = sortrows([lefts; rights])'; % one bup in each col, second row has side bup was on
    #    
    #    adapted = ones(1,size(allbups,2)); 
    #    nclicks = ones(size(adapted)); 
    #    
    #    # now let's go through and figure all the across-side adaptive effects
    #    for c = 1:size(allbups,2)-1,
    #        if allbups(2,c)~=allbups(2,c+1), % if this bup and the next are on opposite sides
    #            dt = allbups(1,c+1) - allbups(1,c);
    #            if dt <= cross_side_suppression,
    #                adapted(c) = 0;
    #                adapted(c+1) = 0;
    #                nclicks(c) = 0.5;
    #                nclicks(c+1) = 0.5;
    #            else
    #                # strength of the cross-side adaptation is weighed by the
    #                # magnitude of the preceeding click
    #                adapted(c+1) = 1 - exp(-dt/tau_psi + log(1 - (adapted(c)*(psi-1)+1)));
    #            end;
    #        end;
    #    end;
    #    
    #    
    #    Lother = real(adapted(allbups(2,:)==-1));
    #    Rother = real(adapted(allbups(2,:)==+1));
    #    Lnclicks = nclicks(allbups(2,:)==-1);
    #    Rnclicks = nclicks(allbups(2,:)==+1);
        println("you didn't implement the cross-stream adaptation system yet!!!!!")
        Lnclicks = ones(size(leftbups));
        Rnclicks = ones(size(rightbups));
    else
        Lnclicks = ones(size(leftbups));
        Rnclicks = ones(size(rightbups));
    end
    
    
    # now take the product of the two effects:
    L = Lsame .* Lother;
    R = Rsame .* Rother;
    return L, R
end


## Adaptation function with both within stream and across stream adaptation
#function make_adapted_cat_clicks(leftbups, rightbups, phi, tau_phi)
#    cross_side_suppression = 0;
#    
#    if abs(phi - 1) > eps()
#        lefts  = [leftbups;  -ones(1,length(leftbups))];
#        rights = [rightbups; +ones(1,length(rightbups))];
#        allbups = sortrows([lefts rights]')'; # one bup in each col, second row has side bup was on
#
#        if length(allbups) <= 1
#            ici = [];
#        else
#            ici = (allbups[1,2:end]  - allbups[1,1:end-1])';
#        end     
#
#        adapted = ones(typeof(phi), 1, size(allbups,2));
#               
#        for i = 2:size(allbups,2)
#            if ici[i-1] <= cross_side_suppression
#                adapted[i-1] = 0;
#                adapted[i] = 0;
#            else
#                last = tau_phi * log(1 - adapted[i-1]*phi);
#                adapted[i] = 1 - exp((-ici[i-1] + last)/tau_phi);
#            end
#        end
#    
#    	adapted = real(adapted);
#    
#    	L = adapted[allbups[2,:] .==-1]';
#    	R = adapted[allbups[2,:] .==+1]';
#    else
#    	# phi was equal to 1, there's no adaptation going on.
#    	L = leftbups;
#    	R = rightbups;
#    end
#
#    return L, R
#end


"""
    make_adapted_cat_clicks(leftbups, rightbups, phi, tau_phi)

    Adaptation function with both within stream and across stream adaptation
    Returns adapted click rates for left and right streams
"""
function make_adapted_cat_clicks(leftbups, rightbups, phi, tau_phi)
    
    if abs(phi - 1) > eps()
        lefts  = [leftbups;  -ones(1,length(leftbups))];
        rights = [rightbups; +ones(1,length(rightbups))];
        if isempty(lefts) & !isempty(rights)
            allbups = rights;
        elseif isempty(rights) & !isempty(lefts)
            allbups = lefts;
        elseif isempty(lefts) & isempty(rights)
            allbups = [];
        else
            allbups = sortslices([lefts rights]',dims=1)'; # one bup in each col, second row has side bup was on
        end

        if length(allbups) <= 1
            ici = [];
            adapted = [];
        else
            ici = (allbups[1,2:end]  - allbups[1,1:end-1])';
        end     

        adapted = ones(typeof(phi), 1, size(allbups,2));
        for i = 2:size(allbups,2)
            if ici[i-1] <= 0
                adapted[i-1] = 0;
                adapted[i] =0;
            else
                #last = tau_phi * log(1 - adapted[i-1]*phi);
                #adapted[i] = 1 - exp((-ici[i-1] + last)/tau_phi);
                adapted[i] = 1+ exp(-ici[i-1]/tau_phi)*(adapted[i-1]*phi -1);
            end
        end
    
    	adapted = real(adapted);
        # Put this in a try catch Loop

        
        if isempty(allbups)
            L = [];
            R = [];
        else
            if any(allbups[2,:] .==-1) 
                L = adapted[allbups[2,:] .==-1]';
            else
                L = [];
            end
            if any(allbups[2,:] .==+1)
                R = adapted[allbups[2,:] .==+1]';
            else
                R = [];
            end
        end
    

    else
    	# phi was equal to 1, there's no adaptation going on.
    	L = leftbups;
    	R = rightbups;
    end
    return L, R
end
