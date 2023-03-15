using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
# Hartree Fock related 
params = Params(w1=96.056,w0=0.7*96.056,vf=2135.4,dθ=1.05π/180)

##
flag = "flavor"
seed = 1
w0s = ["07"]
w0snum = [0.7]
p,q = 1,4
νF = 0+ (1)*p/q
νstr = round(Int,1000*νF)
for w0 in w0s
    metadata = joinpath(fpath,"princeton/data_w$(w0)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    println(load(metadata,"iter_energy")[end])
    println(load(metadata,"iter_err")[end])
    # plot_hf_iterations(metadata)
    plot_spectra(metadata;savename="test.pdf")
end 


## BM basis 
seed = 1
p, q = 1 ,4
flag = "flavor"
νF = 2+ (2)*p/q
νstr = round(Int,1000*νF)
metadata = joinpath(fpath,"princeton/data_w07/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
# plot_density_matrix_bm_valley_spin(metadata)
# plot_density_matrix_bm(metadata)

## strong coupling basis at reference point defined by metadata0

seed = 1
flag  = "flavor"
νF0 = 2+ (0)*p/q
νstr0 = round(Int,1000*νF0)
metadata0 = joinpath(fpath,"princeton/data_w07/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr0).jld2")
# plot_density_matrix_strong_coupling(metadata,metadata0)
plot_density_matrix_strong_coupling_valley_spin(metadata,metadata0)
# 
plot_order_parameters(metadata)
