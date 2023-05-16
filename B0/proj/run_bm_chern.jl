using PyPlot
using JLD2
fpath = pwd()
include(joinpath(fpath,"B0/libs/BMChern_mod.jl"))
include(joinpath(fpath,"B0/libs/plot_helpers.jl"))

phi = 120 #parse(Int,ARGS[1])
strain = 2 #parse(Int,ARGS[2])
# ------------------ Specification ------------------ #
lk = 19
# params = Params(ϵ=0.00,Da=0,dθ=1.06π/180,w1=110,w0=77,vf=2482)
params = Params(ϵ=0.001*strain,Da=-4100,φ=phi*π/180,dθ=1.05π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)
latt = Lattice()
initLattice(latt,params;lk=lk)

bm_path = joinpath(fpath,"feldman/B0/data/strain$(strain)/phi$(phi)/bm_lk$(lk).jld2")

# ------------------ non-interacting part ------------------ #
function compute_bm(latt::Lattice,params::Params;fname::String="placeholder.txt")
    bm = HBM()
    initHBM(bm,latt,params;
            lg=9,_σrotation=false,_calculate_overlap=false,fname=fname)
    return bm
end

bm = compute_bm(latt,params,fname=bm_path);

# ------------------ non-interacting analysis ------------------ #

# kvec = reshape(latt.kvec ./ abs(params.g1),lk,lk)
# tmp = reshape(load(bm_path,"E"),8,8,lk,lk)
# ϵ0 = zeros(Float64,8,lk,lk)
# for i2 in 1:lk, i1 in 1:lk 
#     ϵ0[:,i1,i2] = eigvals(Hermitian(tmp[:,:,i1,i2]))
# end 
# plot_contour_maps(kvec,ϵ0[9,:,:];points=[params.Kt/abs(params.g1)],contourlines=[100.])
# plot_contour_maps(kvec,ϵ0[1,:,:];points=ComplexF64[],contourlines=Float64[],limits=Float64[-10,0])
# iΓ = (lk%2==0) ? (lk÷2) : ((lk-1)÷2+1)
# kcut = real(kvec[:,iΓ])
# Ecut = reshape(ϵ0[1:2:end,:,iΓ],:,length(kcut))
# plot_energy_cuts(kcut,Ecut,lines=Float64[])

# ----------------- plot real space chern states ------------------- #
# check that smooth gauge is properly implemented.
# include(joinpath(fpath,"B0/libs/hybridWannier_mod.jl"));
# hWS_r = initHybridWannierRealSpace(bm); 
# # plot_HybridWannier_RealSpace(hWS_r,params);
# plot_HybridWannier_RealSpaceCombined(hWS_r,params,nk=13);