using PyPlot
using JLD2
fpath = pwd()
include(joinpath(fpath,"B0/libs/HFChern_mod.jl"))
include(joinpath(fpath,"B0/libs/plot_helpers.jl"))

BLAS.set_num_threads(1)

prefix = ARGS[1]
ν = parse(Float64,ARGS[2])
νstr = round(Int,1000*ν)
flag = ARGS[3]
# ------------------ Specification ------------------ #
lk = 23
# params = Params(ϵ=0.00,Da=0,dθ=1.06π/180,w1=110,w0=77,vf=2482)
params = Params(ϵ=0.002,Da=-4100,dθ=1.05π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)
latt = Lattice()
initLattice(latt,params;lk=lk)

bm_path = joinpath(fpath,"feldman/B0/data/bm_lk$(lk).jld2")
hf_path = joinpath(fpath,"feldman/B0/data/$(prefix)_$(flag)_hf_$(νstr)_lk$(lk).jld2")

# ------------------ Hartree-Fock part ------------------ #
function compute_hf(ν::Float64,latt::Lattice,params::Params;fname::String="placeholder.txt")
    hf = HartreeFock()
    run_HartreeFock(hf,params,latt,bm_path; ν=ν,_Init=flag,savename=fname)
    return hf
end

hf = compute_hf(ν,latt,params,fname=hf_path);