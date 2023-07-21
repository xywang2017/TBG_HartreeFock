using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))

#
# Info and folder name
# ------------------------------------------------------------------------------ # 
twist_angle = 124
foldername = "$(twist_angle)_strain"
params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)

# ----------------------------------Hartree Fock spectrum-------------------------------------------- # 
flag = "random"
seed = 3
w0s = ["07"]
w0snum = [0.7]
p,q = 1,4
νF = (-3)+(-1)*p/q
νstr = round(Int,1000*νF)
hf = 0
for w0 in w0s
    metadata = joinpath(fpath,"$(foldername)/B/data_w$(w0)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    println("HF energy: ",load(metadata,"iter_energy")[end])
    println("Convergence: ",load(metadata,"iter_err")[end])
    plot_spectra(metadata;savename="test.pdf")
end 

# -----------------------------------Density matrix analysis ------------------------------------------- # 
seed = 1
p, q = 1,4
flag = "bm_cascade"
νF = (-3) + (-1)*p/q
νstr = round(Int,1000*νF)
metadata = joinpath(fpath,"124_strain/B/data_w07/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
# plot_density_matrix_bm_valley_spin(metadata)
plot_density_matrix_bm(metadata,ik=2)
# plot_density_matrix_sublattice(metadata)
plot_density_matrix_sublattice_full(metadata)


## strong coupling basis at reference point defined by metadata0
seed = 1
flag  = "random"
νF0 = 0+ (0)*p/q
νstr0 = round(Int,1000*νF0)
metadata0 = joinpath(fpath,"princeton/data_w07/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr0).jld2")
# plot_density_matrix_strong_coupling(metadata,metadata0)
# plot_density_matrix_strong_coupling_valley_spin(metadata,metadata0)
# 
# plot_order_parameters(metadata)