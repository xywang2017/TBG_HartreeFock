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
νF = 0+ (2)*p/q
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
flag = "random"
νF = 2+ (2)*p/q
νstr = round(Int,1000*νF)
metadata = joinpath(fpath,"princeton/data_w07/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
# plot_density_matrix_bm_valley_spin(metadata)
# plot_density_matrix_bm(metadata)
# plot_density_matrix_sublattice(metadata)

## strong coupling basis at reference point defined by metadata0

seed = 1
flag  = "flavor"
νF0 = 0+ (0)*p/q
νstr0 = round(Int,1000*νF0)
metadata0 = joinpath(fpath,"princeton/data_w07/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr0).jld2")
# plot_density_matrix_strong_coupling(metadata,metadata0)
plot_density_matrix_strong_coupling_valley_spin(metadata,metadata0)
# 
plot_order_parameters(metadata)


##
x = [0,1,2,3,4,5,6,7]
y = [13.23;5.05;7.36;8.26;5.92;4.72;16.44;1.92];
fig = figure(figsize=(4,3))
plot(x,y,"b^")
xticks(x,["(0,0)","(0,1)","(0,2)","(0,3)","(0,4)","(2,0)","(2,1)","(2,2)"])
ylim([0,19])
yticks(0:4:20)
xlabel("(s,t)")
ylabel("Δ (meV)")
tight_layout()
savefig("gaps_strongcoupling.pdf")
display(fig)
close(fig)