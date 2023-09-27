using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))

#
# Info and folder name
# ------------------------------------------------------------------------------ # 
twist_angle = 132
foldername = "zeeman/$(twist_angle)_strain"
params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)

# ----------------------------------Hartree Fock spectrum-------------------------------------------- # 
s,t = -3,-1
p,q = 3,11
νF = (s)+(t)*p/q
println(νF)
νstr = round(Int,1000*νF)
metadata = joinpath(fpath,"$(foldername)/_$(p)_$(q)/1_random_tL_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
if !isfile(metadata)
    metadata = joinpath(fpath,"$(foldername)/_$(p)_$(q)/1_bm_cascade_tL_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
end
if isfile(metadata)
    E = load(metadata,"iter_energy")[end]
    for flag in ["flavor","random","random_tL","bm","strong","bm_cascade","bm_cascade_tL"], seed in 1:12
        metadata0 = joinpath(fpath,"$(foldername)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
        if isfile(metadata0)
            E0 = load(metadata0,"iter_energy")[end]
            if E0<=E 
                E, metadata = E0, metadata0 
            end
        end
    end
    println(metadata)
end
println("HF energy: ",load(metadata,"iter_energy")[end])
println("Convergence: ",load(metadata,"iter_err")[end])

# -----------------------------------Density matrix analysis ------------------------------------------- # 
# plot_spectra(metadata;savename="test.png")
# plot_density_matrix_bm_valley_spinv0(metadata,ik=1,savename="tmp3.png",jj=3)
# plot_density_matrix_bm(metadata,ik=1)
plot_density_matrix_bm_half(metadata,ik=1)
plot_density_matrix_global_order_parameters(metadata)

# plot_density_matrix_bm(metadata,ik=1)
# plot_density_matrix_sublattice(metadata)
# plot_density_matrix_sublattice_full(metadata)

# plot_density_matrix_global_order_parameters(metadata)
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
