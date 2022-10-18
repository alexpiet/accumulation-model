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

    # load rat data
    println(res_dir)
    println(ratname)
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
    nt = length(data["pokedR"])
    println(string(nt)*" trials in this dataset")
    # evaluate model LL just to make sure its correct
    NLL = compute_LL(data, params);
    println(NLL)
    bad_NLL = abs(NLL -  fit["f"]) > 1
    if bad_NLL
        println("Oh no! Julia version is not within tolerance of MATLAB values")
        temp = NLL - fit["f"];
        println(temp)
    else
        println("Good. Julia NLL is within tolerance of MATLAB values")
    end

    # compute hessian using autodiff
    autodiff_hessian = ForwardDiff.hessian(x->compute_LL(data,x), params)

    # save new hessian
    matwrite(save_path,Dict("autodiff_hessian"=>autodiff_hessian,"julia_nll"=>NLL,"bad_nll"=>bad_NLL))
    println("saved data")

    return autodiff_hessian

end





