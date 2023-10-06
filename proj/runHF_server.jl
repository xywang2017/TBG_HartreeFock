using JLD2
using Printf
using LinearAlgebra
fpath = pwd()
# include(joinpath(fpath,"libs/MagneticFieldHF_tLSymmetric.jl"))
# include(joinpath(fpath,"libs/MagneticFieldHF.jl"))

BLAS.set_num_threads(1)
# dir =  "/Volumes/Data/Code/TBG_HartreeFock/"
dir = ""
#
## Hartree Fock related 
p = parse(Int,ARGS[1])
q = parse(Int,ARGS[2])
ν = parse(Float64,ARGS[3])
νstr = round(Int,1000*ν)
flag = ARGS[4]
w0 = ARGS[5]
seed = ARGS[6]
twist_angle = parse(Float64,ARGS[7])
_is_strain = ARGS[8]
_is_symmetric = ARGS[9]

if isequal(_is_symmetric,"symmetric")
    include(joinpath(fpath,"libs/MagneticFieldHF_tLSymmetric.jl"))
else
    include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
end

foldername = @sprintf "%d_%s" round(Int,twist_angle*100) _is_strain
if ! isdir(joinpath(fpath,"zeeman/$(foldername)/_$(p)_$(q)"))
    mkpath(joinpath(fpath,"zeeman/$(foldername)/_$(p)_$(q)"))
end

if isequal(_is_symmetric,"symmetric")
    savename = joinpath(fpath,"zeeman/$(foldername)/_$(p)_$(q)/$(seed)_$(flag)_tL_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
else
    savename = joinpath(fpath,"zeeman/$(foldername)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
end

ϕ = p//q 
if isequal(flag,"flavor")
    _Init = "Flavor U(4)"
elseif isequal(flag,"chern")
    _Init = "Sublattice"
elseif isequal(flag,"random")
    _Init = "Random"
elseif isequal(flag,"bm")
    _Init = "bm"
elseif isequal(flag,"bm_cascade")
    _Init = "bm_cascade"
elseif isequal(flag,"vssymmetric")
    _Init = "vssymmetric"
elseif isequal(flag,"strong")
    _Init = " "
end

println("Running parameters: ","ϕ=",ϕ,", ν=",ν,", Init=",flag,", w0=",w0)
println(savename)

if isequal(_is_strain,"strain")
    params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*π/180,w1=110,w0=parse(Float64,ARGS[5])*0.1*110,vf=2482)
elseif isequal(_is_strain,"nostrain")
    params = Params(ϵ=0.00,Da=0.0,φ=0.0*π/180,dθ=twist_angle*π/180,w1=110,w0=parse(Float64,ARGS[5])*0.1*110,vf=2482)
end
initParamsWithStrain(params)
hf = HartreeFock()

if !isequal(flag,"strong")
    iter_err, iter_energy = run_HartreeFock(hf,params,ν=ν,ϕ=ϕ,prefix=dir*"NonInt/$(foldername)/",_Init=_Init,savename=savename)
else
    savename0 = dir*"zeeman/112_$(_is_strain)/_$(p)_$(q)/1_random_tL_init_HF_1_8_nu_-1375.jld2"
    hf0 = load(savename0,"hf")
    P0,H0 = hf0.P,hf0.H
    iter_err, iter_energy = run_HartreeFock(hf,params,ν=ν,ϕ=ϕ,prefix=dir*"NonInt/$(foldername)/",_Init=" ",H0=H0,P0=P0,savename=savename)
end
