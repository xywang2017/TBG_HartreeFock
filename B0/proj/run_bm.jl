using PyPlot
using JLD2
fpath = pwd()
include(joinpath(fpath,"B0/libs/BM_mod.jl"))
include(joinpath(fpath,"B0/libs/plot_helpers.jl"))

# ------------------ Specification ------------------ #
lk = 65
# params = Params(ϵ=0.002,Da=-4100,dθ=1.06π/180,w1=110,w0=77,vf=2482)
# params = Params(ϵ=0.002,Da=-4100,dθ=1.27π/180,w1=110,w0=77,vf=2680)
params = Params(ϵ=0.003,Da=-4100,dθ=1.05π/180,w1=110,w0=77,vf=2482)
# params = Params(ϵ=0.00,Da=0,w1=96.056,w0=0.7*96.056,vf=2135.4,dθ=1.05π/180)
initParamsWithStrain(params)
latt = Lattice()
initLattice(latt,params;lk=lk)

bm_path = joinpath(fpath,"feldman/B0/data/test_bm_lk$(lk).jld2")

# ------------------ non-interacting part ------------------ #
function compute_bm(latt::Lattice,params::Params;fname::String="placeholder.txt")
    bm = HBM()
    initHBM(bm,latt,params;
            lg=9,_σrotation=true,_calculate_overlap=false,fname=fname)
    return bm
end

bm = compute_bm(latt,params,fname=bm_path);

# ------------------ non-interacting analysis ------------------ #

kvec = reshape(latt.kvec ./ abs(params.g1),lk,lk)
ϵ0 = reshape(load(bm_path,"E"),:,lk,lk)
# plot_contour_maps(kvec,ϵ0[9,:,:];points=[params.Kt/abs(params.g1)],contourlines=[100.])
plot_contour_maps(kvec,ϵ0[9,:,:];points=ComplexF64[],contourlines=Float64[])
iΓ = (lk%2==0) ? (lk÷2) : ((lk-1)÷2+1)
kcut = real(kvec[:,iΓ])
Ecut = reshape(ϵ0[1:2:end,:,iΓ],:,length(kcut))
plot_energy_cuts(kcut,Ecut,lines=Float64[])