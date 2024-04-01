"""
    compute_hessian(ratname; res_dir="./", overwrite=true)

    Compute the hessian of the likelihood function relation for the choice model predicting right 
    and left pokes from clicks times and model parameters at the maximum likelihood estimate for a given rat

"""
function compute_hessian(ratname; res_dir="./", overwrite=true)

    # Loop over rats if provided as group
    if ratname isa Array
        nrats   = length(ratname)
        hess    = Array{Any, 1}(undef,nrats)
        for rr  = 1:nrats
            hess[rr] = compute_hessian(ratname[rr], res_dir=res_dir, overwrite=overwrite)
        end
        return hess
    end

    println("\nComputing hessian for "*ratname)
    println("----------------------------------")
    
    # load rat data
    fitfn = joinpath(res_dir,"fit_analytical_"*ratname*".mat");
    save_path = joinpath(res_dir,"julia_hessian_"*ratname*".mat")

    # load hess if it exists
    if isfile(save_path) & !overwrite
        println("already saved hessian for "*ratname) 
        hess = matread(save_path)
        return hess
    end

    # load fit if it exists
    if isfile(fitfn)
        fit_data = matread(fitfn);
        fit      = fit_data["fit"] 
    else
        error("couldn't find fit file "*fitfn)
    end
    println("Analytical choice model fit loaded")
    fit     = fit_data["fit"];
    println("Contains field 'fit'")
    data    = fit_data["data"]
    pokedR = data["avgdata"]["pokedR"]
    if haskey(data,"rawdata")
        data = data["rawdata"]
    end
    data["pokedR"]=(pokedR.==1)
    nt = length(pokedR)

    println("Contains field 'data' with "*string(nt)*" trials")
    params  = fit["final"];
    println("Contains field 'final' with parameters and priors")

    # Figure out whether to set priors on model parameters
    prior_mean = fit["prior_mean"]
    prior_var = fit["prior_var"]

    # Print the parameters and priors that will be used
    for i = 1:length(params)
        # print the parameters to 3 decimals
        println(fit["param_names"][i]*": "*string(round(params[i],digits=3))*" (prior mean: "*string(round(prior_mean[i],digits=3))*")")
    end

    # Compute model likelihood for these parameters and priors
    NLL  = compute_LL(data, params; prior_var=prior_var, prior_mean=prior_mean);

    # Report the NLL 
    println("MATLAB NLL is: "*string(round(fit["f"],digits=3))*" Julia NLL is: "*string(round(NLL,digits=3)))

    # Check if NLL is within tolerance of MATLAB values
    bad_NLL = abs(NLL -  fit["f"]) > 1
    if bad_NLL
        println("Oh no! Julia version is not within tolerance of MATLAB values")
        temp = NLL - fit["f"];
        println(temp)
    else
        println("Good. Julia NLL is within tolerance of MATLAB values. (Julia - MATLAB = "*string(NLL - fit["f"])*")")
    end

    # compute hessian using autodiff
    autodiff_hessian = ForwardDiff.hessian(x->compute_LL(data, x; prior_mean = prior_mean, prior_var=prior_var), params)
    
    res = Dict("NLL"=>NLL, "autodiff_hessian"=>autodiff_hessian, "bad_NLL"=>bad_NLL)

    return res

end


function save_hessian(ratname; res_dir="./", overwrite=true)

    # Loop over rats if provided as group
    if ratname isa Array
        nrats   = length(ratname)
        res    = Array{Any, 1}(undef,nrats)
        for rr  = 1:nrats
            res[rr] = save_hessian(ratname[rr], res_dir=res_dir, overwrite=overwrite)
        end
        return res
    end

    res = compute_hessian(ratname, res_dir=res_dir, overwrite=overwrite)
    save_path = joinpath(res_dir,"julia_hessian_"*ratname*".mat")
    matwrite(save_path, res)
    println("saved data in "*save_path)

    return res
end


