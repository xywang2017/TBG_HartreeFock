using PyPlot
using JLD2
fpath = pwd()
include(joinpath(fpath,"B0/libs/BM_mod.jl"))
include(joinpath(fpath,"B0/libs/plot_helpers.jl"))

# ------------------ Specification ------------------ #
lk = 33
params = Params(ϵ=0.00,φ=0.0,Da=-4100,dθ=1.38π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)
latt = Lattice()
initLattice(latt,params;lk=lk)

bm_path = joinpath(fpath,"B0/data/bm_lk$(lk).jld2")

# ------------------ non-interacting part ------------------ #
function compute_bm(latt::Lattice,params::Params;fname::String="placeholder.txt")
    bm = HBM()
    initHBM(bm,latt,params;
            lg=9,_σrotation=false,_calculate_overlap=false,fname=fname)
    return bm
end

bm = compute_bm(latt,params,fname=bm_path);

# layer_P= real( calculateLayerPolarization(bm) );

# ------------------ non-interacting analysis ------------------ #

kvec = reshape(latt.kvec ./ abs(params.g1),lk+1,lk+1)
# kvec = reshape(latt.k1,:,1) .+ 1im*reshape(latt.k2,1,:) 
ϵ0 = reshape(load(bm_path,"E"),:,lk+1,lk+1)
# layer = reshape(layer_P,2,2,lk,lk)
# plot_contour_maps(kvec,ϵ0[9,:,:];points=[params.Kt/abs(params.g1)],contourlines=[100.])
plot_contour_maps(kvec,ϵ0[5,:,:];points=[params.Kt/abs(params.g1)],contourlines=Float64[])
# plot_contour_maps(kvec,layer[2,2,:,:];points=kvec[k_indx],contourlines=Float64[])
iΓ = (lk%2==0) ? (lk÷2) : ((lk-1)÷2+1)
kcut = real(kvec[:,iΓ])
Ecut = reshape(ϵ0[1:2:end,:,iΓ],:,length(kcut))
plot_energy_cuts(kcut,Ecut,lines=Float64[])

minimum(ϵ0[13,:,:]) - maximum(ϵ0[12,:,:])
maximum(ϵ0[12,:,:])

ee = 1.6e-19
ϵϵ = 8.8541878128e−12	
aa = 2.46e-10
ϵr = 15.0
Vcoulomb = ee/(4π*ϵϵ*ϵr* sqrt(abs(params.a1)*abs(params.a1))*aa) * 1e3 * exp(-1/4)