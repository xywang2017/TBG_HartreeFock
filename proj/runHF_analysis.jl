using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))

#
# Info and folder name
# ------------------------------------------------------------------------------ # 
twist_angle = 120
foldername = "$(twist_angle)_strain"
params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)

# ----------------------------------Hartree Fock spectrum-------------------------------------------- # 
# flag = "random"
# seed = 2
# for sts in [[0,-4],[-1,-3],[-2,-2],[-3,-1]]
    # s,t = sts[1], sts[2]
    s,t = -1,-3
    p,q = 1,4
    νF = (s)+(t)*p/q
    νstr = round(Int,1000*νF)
    metadata = joinpath(fpath,"$(foldername)/_$(p)_$(q)/1_random_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    if !isfile(metadata)
        metadata = joinpath(fpath,"$(foldername)/_$(p)_$(q)/1_bm_cascade_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    end
    if isfile(metadata)
        E = load(metadata,"iter_energy")[end]
        for flag in ["flavor","random","chern","bm","strong","bm_cascade"], seed in 1:10
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
    # plot_spectra(metadata;savename="test.pdf")
    plot_density_matrix_bm_valley_spinv2(metadata,ik=12,savename="124_DensityMat_HFM_$(s)_$(t).png")
# end
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
