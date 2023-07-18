using PyPlot
using JLD2
using Printf
fpath = pwd()
include(joinpath(fpath,"libs/HFChern_mod.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))

BLAS.set_num_threads(1)

prefix = ARGS[1]
ν = parse(Float64,ARGS[2])
νstr = round(Int,1000*ν)
flag = ARGS[3]
twist_angle = parse(Float64,ARGS[4])
_is_strain = ARGS[5]
lk = parse(Int,ARGS[6])
foldername =@sprintf "%d_%s" round(Int,twist_angle*100) _is_strain 
# ------------------ Specification ------------------ #
if isequal(_is_strain,"strain")
    params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*π/180,w1=110,w0=77,vf=2482)
elseif isequal(_is+strain,"nostrain")
    params = Params(ϵ=0.00,Da=0.0,φ=0.0*π/180,dθ=twist_angle*π/180,w1=110,w0=77,vf=2482)
end
initParamsWithStrain(params)
latt = Lattice()
initLattice(latt,params;lk=lk)

bm_path = joinpath(fpath,"$(foldername)/B0/bm_lk$(lk).jld2")
hf_path = joinpath(fpath,"$(foldername)/B0/$(prefix)_$(flag)_hf_$(νstr)_lk$(lk).jld2")

println(hf_path)
# ------------------ Hartree-Fock part ------------------ #
function compute_hf(ν::Float64,latt::Lattice,params::Params;fname::String="placeholder.txt")
    hf = HartreeFock()
    run_HartreeFock(hf,params,latt,bm_path; ν=ν,_Init=flag,savename=fname)
    return hf
end

hf = compute_hf(ν,latt,params,fname=hf_path);
