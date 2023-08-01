using PyPlot
using Printf
using JLD2
fpath = pwd()
include(joinpath(fpath,"B0/libs/BMChern_mod.jl"))
include(joinpath(fpath,"B0/libs/plot_helpers.jl"))

twist_angle = 1.28 #parse(Float64,ARGS[1])
_is_strain = "strain" #ARGS[2]
lk = 15 #parse(Int,ARGS[3])

foldername =@sprintf "%d_%s" round(Int,twist_angle*100) _is_strain 
# ------------------ Specification ------------------ #
if isequal(_is_strain,"strain")
    params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*π/180,w1=110,w0=77,vf=2482)
elseif isequal(_is_strain,"nostrain")
    params = Params(ϵ=0.000,Da=0.0,φ=0.0*π/180,dθ=twist_angle*π/180,w1=110,w0=77,vf=2482)
end
initParamsWithStrain(params)
latt = Lattice()
initLattice(latt,params;lk=lk)

bm_path = joinpath(fpath,"$(foldername)/B0/bm_lk$(lk).jld2")

# if !isdir(bm_path)
#     mkpath(bm_path)
# end

# ------------------ non-interacting part ------------------ #
function compute_bm(latt::Lattice,params::Params;fname::String="placeholder.txt")
    bm = HBM()
    initHBM(bm,latt,params;
            lg=9,_σrotation=false,_calculate_overlap=true,fname=fname)
    return bm
end

bm = compute_bm(latt,params,fname=bm_path);
