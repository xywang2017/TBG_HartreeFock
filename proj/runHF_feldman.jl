using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
## Hartree Fock related 
params = Params()
initParamsWithStrain(params)

##
flag = "random"
seed = 2
w0cs = Float64[]
w0s = ["07"]
w0snum = [0.7]
σz = []
p,q = 1, 6
νF = 2+ (2)*p/q
νstr = round(Int,1000*νF)
for w0 in w0s
    metadata = joinpath(fpath,"feldman/data_w$(w0)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    hf = load(metadata,"hf")
    # metadata = "typical_starting_point.jld2"
    println(load(metadata,"iter_energy")[end])
    println(load(metadata,"iter_err")[end])
    # plot_hf_iterations(metadata)
    # ϵk = load(metadata,"spectrum")
    # σzτz = load(metadata,"chern")
    actual_ν = 8*round(Int,(νF+4)/8*length(hf.ϵk)) / length(hf.ϵk)-4
    ϵk = hf.ϵk 
    σzτz = hf.σzτz
    push!(σz, plot_spectra(ϵk,σzτz,actual_ν,params;savename="test.pdf"))
end 

## BM basis 
seed = 1
νF = 1+ (3)*p/q
νstr = round(Int,1000*νF)
metadata = joinpath(fpath,"feldman/data_w07/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
# plot_density_matrix_bm_valley_spin(metadata)
plot_density_matrix_bm(metadata)

## strong coupling basis at reference point defined by metadata0

νF0 = 0+ (0)*p/q
νstr0 = round(Int,1000*νF0)
metadata0 = joinpath(fpath,"feldman/data_w07/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr0).jld2")
plot_density_matrix_strong_coupling(metadata,metadata0)
plot_density_matrix_strong_coupling_valley_spin(metadata,metadata0)
# 
plot_order_parameters(metadata)
