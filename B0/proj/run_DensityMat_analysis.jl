using PyPlot
using JLD2
fpath = pwd() 
include(joinpath(fpath,"B0/libs/DensityMatrix_reduction.jl"))
include(joinpath(fpath,"B0/libs/plot_helpers.jl"))

# ------------------ Load Hartree Fock results ------------------ #
prefix =1
flag = "random"
twist_angle = 1.20
_is_strain = "strain"
foldername = @sprintf "%d_%s" round(Int,twist_angle*100) _is_strain
ν = -2.0
νstr = round(Int,1000*ν)
lk = 13
params = Params(ϵ=0.002,Da=-4100,dθ=twist_angle*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)
latt = Lattice()
initLattice(latt,params;lk=lk)
kvec = reshape(latt.kvec ./ abs(params.g1),lk,lk)
# kvec = reshape(latt.k1,:,1) .+ 1im*reshape(latt.k2,1,:)

hf_path = joinpath(fpath,"$(foldername)/B0/$(prefix)_$(flag)_hf_$(νstr)_lk$(lk).jld2")
hf = load(hf_path,"hf");
println("HF energy: ",load(hf_path,"iter_energy")[end])
println("HF convergence: ",load(hf_path,"iter_err")[end])

# ---------------- Density Matrix analysis ---------------- # 
dm = constructDensityMat(hf);
# checkReconstructionValidity(dm,collect(60:64))
plot_contour_maps(kvec,reshape(dm.δs[4,:],lk,lk),points=ComplexF64[0+0im],contourlines=[100.],limits=Float64[])

# ----------------------- plot correlation values --------  #
# plot_corr_values(dm)
# ------------------- Analysis of structures of 8x8 matrices Oϕs ----------------- # 
for idx in 60:64
    # plot_formfactor_info(dm,idx)
    plot_formfactor_info_band_basis(dm,idx)
end
