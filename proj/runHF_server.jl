using JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))

BLAS.set_num_threads(1)
#
## Hartree Fock related 
p = parse(Int,ARGS[1])
q = parse(Int,ARGS[2])
ν = parse(Float64,ARGS[3])
νstr = round(Int,1000*ν)
flag = ARGS[4]
w0 = ARGS[5]
seed = ARGS[6]
savename = joinpath(fpath,"feldman/B/data_w$(w0)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
ϕ = p//q 
if isequal(flag,"flavor")
    _Init = "Flavor U(4)"
elseif isequal(flag,"chern")
    _Init = "Sublattice"
elseif isequal(flag,"random")
    _Init = "Random"
elseif isequal(flag,"bm")
    _Init = "bm"
elseif isequal(flag,"strong")
    _Init = " "
end

println("Running parameters: ","ϕ=",ϕ,", ν=",ν,", Init=",flag,", w0=",w0)
println(savename)
# params = Params(w1=96.056,w0=parse(Float64,ARGS[5])*0.1*96.056,vf=2135.4,dθ=1.05π/180)
params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=1.05π/180,w1=110,w0=parse(Float64,ARGS[5])*0.1*110,vf=2482)
# params = Params(ϵ=0.00,Da=-4100,φ=0.0*π/180,dθ=1.2π/180,w1=110,w0=parse(Float64,ARGS[5])*0.1*110,vf=2482)
initParamsWithStrain(params)
hf = HartreeFock()

if !isequal(flag,"strong")
    iter_err, iter_energy = run_HartreeFock(hf,params,ν=ν,ϕ=ϕ,prefix="feldman/B/data_w$(w0)/_$(p)_$(q)/",_Init=_Init,savename=savename)
else
    savename0 = joinpath(fpath,"feldman/B/data_w$(w0)/_$(p)_$(q)/$(seed)_random_init_HF_$(p)_$(q)_nu_0.jld2")
    hf0 = load(savename0,"hf")
    P0,H0 = hf0.P,hf0.H
    iter_err, iter_energy = run_HartreeFock(hf,params,ν=ν,ϕ=ϕ,prefix="feldman/B/data_w$(w0)/_$(p)_$(q)/",_Init=" ",H0=H0,P0=P0,savename=savename)
end

# save(savename,"H",hf.H,"P",hf.P,"spectrum",hf.ϵk,"chern",hf.σzτz,"iter_err",iter_err,"iter_energy",iter_energy)
#