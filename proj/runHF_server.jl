using JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
#
## Hartree Fock related 
p = parse(Int,ARGS[1])
q = parse(Int,ARGS[2])
ν = parse(Float64,ARGS[3])
νstr = round(Int,1000*ν)
flag = ARGS[4]
w0 = ARGS[5]
seed = ARGS[6]
savename = joinpath(fpath,"princeton/data_w$(w0)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
ϕ = p//q 
if isequal(flag,"flavor")
    _Init = "Flavor U(4)"
elseif isequal(flag,"chern")
    _Init = "Sublattice"
elseif isequal(flag,"random")
    _Init = "Random"
elseif isequal(flag,"strong")
    _Init = " "
end

println("Running parameters: ","ϕ=",ϕ,", ν=",ν,", Init=",flag,", w0=",w0)
println(savename)
params = Params(w1=96.056,w0=parse(Float64,ARGS[5])*0.1*96.056,vf=2135.4,dθ=1.05π/180)
# params = Params(w1=110.0,w0=parse(Float64,ARGS[5])*0.1*110)
# initParamsWithStrain(params)
hf = HartreeFock()
iter_err, iter_energy = run_HartreeFock(hf,params,ν=ν,ϕ=ϕ,prefix="princeton/data_w$(w0)/_$(p)_$(q)/",_Init=_Init,savename=savename)

# savename0 = joinpath(fpath,"princeton/data_w$(w0)/_$(p)_$(q)/$(seed)_flavor_init_HF_$(p)_$(q)_nu_0.jld2")
# hf0 = load(savename0,"hf")
# P0,H0 = hf0.P,hf0.H0
# iter_err, iter_energy = run_HartreeFock(hf,params,ν=ν,ϕ=ϕ,prefix="princeton/data_w$(w0)/_$(p)_$(q)/",_Init=" ",H0=H0,P0=P0)

# save(savename,"H",hf.H,"P",hf.P,"spectrum",hf.ϵk,"chern",hf.σzτz,"iter_err",iter_err,"iter_energy",iter_energy)
#