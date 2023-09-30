using PyPlot,JLD2
using Printf 
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))

#
# Info and folder name
# ------------------------------------------------------------------------------ # 
twist_angles = [105; collect(106:2:138)] 
twist_angle = 105
# for twist_angle in twist_anglest
foldername = "zeeman/$(twist_angle)_strain"
params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)

# ----------------------------------Hartree Fock spectrum-------------------------------------------- # 
s,t = 0,-4
p,q = 2,5
νF = (s)+(t)*p/q
νstr = round(Int,1000*νF)
metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="nu_$(νstr)",_printinfo=true)

# -----------------------------------Density matrix analysis ------------------------------------------- # 
# plot_spectra(metadata;savename="test.png")
str = @sprintf "θ=%.2f" 0.01*twist_angle
# plot_density_matrix_bm_valley_spinv0(metadata,ik=1,savename="$(str).png",jj=3,titlestr=str)
# plot_density_matrix_sublattice_v0(metadata,titlestr=str)
# end
plot_density_matrix_bm(metadata,ik=1)
# plot_density_matrix_bm_half(metadata,ik=1)
# test_tL2_breaking(metadata)
# plot_density_matrix_global_order_parameters(metadata)

# plot_density_matrix_sublattice(metadata)
# plot_density_matrix_sublattice_full(metadata)
# plot_density_matrix_valley_spin_density_tL2(metadata)

## strong coupling basis at reference point defined by metadata0
seed = 1
flag  = "random"
νF0 = 0+ (0)*p/q
νstr0 = round(Int,1000*νF0)
metadata0 = joinpath(fpath,"$(foldername)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr0).jld2")
# plot_density_matrix_strong_coupling(metadata,metadata0)
# plot_density_matrix_strong_coupling_valley_spin(metadata,metadata0)
# 
# plot_order_parameters(metadata)
