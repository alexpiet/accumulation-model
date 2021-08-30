# a script for computing the model Hessian at the MLE point using autodiff
# does not actually perform optimization

# set up modules for analysis
using MAT
using ForwardDiff
using SpecialFunctions
using Glob
include("/Users/oroville/projects/pbups_dyn/code/accumulation-model/compute_LL.jl")

# move to directory with the fits
res_dir = "/users/oroville/projects/pbups_dyn/results/"
res_dir = "/users/oroville/projects/waiting/results/"
#cd(res_dir)
# which model fitting group to look at
group = "ephys";
rats = ["H037", "H066", "H084", "H129", "H140", "H191"];
rats = ["Z255", "Z256", "Z258", "J244"];
# which rats to use, all with H0- prefix
##rats = [33 37 38 39 40 43 45 58 61 65 66 67 83 84];
##rats = ["B052", "B053", "B065", "B069", "B074", "B083", "B090", "B093", "B097", "B102", "B103", "B104", "B105", "B106", "B107", "B111", "B112", "B113", "B115"];
##rats = ["H034b", "H034d","H036b","H036d","H036a","H045b","H045d","H045a","H046b","H046d","H046a"];
##rats = ["H034d1","H034d2", "H036d1","H036d2","H036a1", "H036a2", "H045d1","H045d2","H045a1", "H045a2", "H046d1","H046d2","H046a1", "H046a2", "metaa1","metaa2","metaa", "metad1","metad2","metad", "metab"];
##rats = ["metad2_c1","metad2_c2","metad2_c3"];
#rats = ["H065","H067"];
#rats  = ["H191" "H176"]


#filestr = [ "H153_tofit_40Hz" "H153_tofit_20and40Hz"];
#filestr = glob("fit_ana*")
#filestr = ["fit_analytical_H126_tofit_8and40Hz.mat"]

#filestr = ["Z255", "J244"]
#filestr = ["Z255"]

filestr = rats

overwrite = true;
# Iterate over rats
for i = 1:length(filestr) 
    # load model parameters at MLE point from matlab
    file_suffix = filestr[i]
    fn = res_dir*"fit_analytical_"*file_suffix*".mat";
    save_path = res_dir*"julia_hessian_"*file_suffix*".mat"
    println(fn)
    if isfile(save_path) & !overwrite
        println("already saved hessian"*file_suffix) 
        continue
    end
    println("Now analyzing rat "*file_suffix)
    if ~isfile(fn)
        error("couldn't find that file")
    end
    #try
        println(fn)
        fit_data= matread(fn);
        println("fit loaded")
        fit     = fit_data["fit"];
        println("found field fit")
        data    = fit_data["data"]
        nt = length(data["pokedR"])
        println(string(nt)*" trials in this dataset")
        println("found field data")
        params  = fit["final"];
        println("found field final")
    #catch
    #    println("problem loading file")
    #end
    nt = length(data["pokedR"])
    println(string(nt)*" trials in this dataset")
    # evaluate model LL just to make sure its correct
    NLL = compute_LL(data, params);
    println(NLL)
    if abs(NLL -  fit["f"]) > 1
        println("Oh no! Julia version is not within tolerance of MATLAB values")
        temp = NLL - fit["f"];
        println(temp)
    else
        println("Julia version with tolerance of MATLAB values")
    end

    # compute hessian using autodiff
    autodiff_hessian = ForwardDiff.hessian(x->compute_LL(data,x), params)

    # save new hessian
    matwrite(save_path,Dict([("autodiff_hessian",autodiff_hessian)]))
    println("saved data")
   # matwrite(string("julia_hessian_", $file_suffix ".mat"),Dict([("autodiff_hessian",autodiff_hessian)]))
end



