using JLD2
fpath = pwd()
# include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/MagneticFieldHFv1.jl"))
#
## Hartree Fock related 
p = parse(Int,ARGS[1])
q = parse(Int,ARGS[2])
ν = parse(Float64,ARGS[3])
νstr = round(Int,1000*ν)
flag = ARGS[4]
w0 = ARGS[5]
savename = joinpath(fpath,"feldman/data_w$(w0)/1_$(p)_$(q)/_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
ϕ = p//q 
if isequal(flag,"flavor")
    _Init = "Flavor U(4)"
elseif isequal(flag,"chern")
    _Init = "Sublattice"
elseif isequal(flag,"random")
    _Init = "Random"
end

println("Running parameters: ","ϕ=",ϕ,", ν=",ν,", Init=",flag,", w0=",w0)
println(savename)
# params = Params(w1=96.056,w0=parse(Float64,ARGS[5])*0.1*96.056)
params = Params()
initParamsWithStrain(params)
hf = HartreeFock()
iter_err, iter_energy = run_HartreeFock(hf,params,ν=ν,ϕ=ϕ,prefix="feldman/data_w$(w0)/_$(p)_$(q)/",_Init=_Init,savename=savename)

# P0 = load(savename,"P")
# iter_err, iter_energy = run_HartreeFock(hf,params,ν=ν,ϕ=ϕ,prefix="feldman/data_w$(w0)/_$(p)_$(q)/",_Init=" ",P0=P0)

# save(savename,"H",hf.H,"P",hf.P,"spectrum",hf.ϵk,"chern",hf.σzτz,"iter_err",iter_err,"iter_energy",iter_energy)
#