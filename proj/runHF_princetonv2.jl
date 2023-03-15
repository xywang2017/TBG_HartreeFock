using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
# Hartree Fock related 
params = Params()
initParamsWithStrain(params)


##
flag = "flavor"
seed = 1
w0s = ["07"]
w0snum = [0.7]
p,q = 1,4
νF = 2 + (1)*p/q
νstr = round(Int,1000*νF)
for w0 in w0s
    metadata = joinpath(fpath,"princeton/data_w$(w0)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    println(load(metadata,"iter_energy")[end])
    println(load(metadata,"iter_err")[end])
    # plot_hf_iterations(metadata)
    plot_spectra(metadata;savename="test.pdf")
end 

# ------------------ plot_spectra_collectively_at_different_flux --------------- # 
ϕs = [1//2;5//12;2//5;1//3;2//7;1//4;1//5;1//6;1//8;1//10;1//12;1//14]
w0 = "07"
metadatas = String[]
for ϕ in ϕs 
    flag, seed = "random", 5
    p,q = numerator(ϕ), denominator(ϕ)
    νstr = round(Int,1000*(1+3*p/q))
    if q==12 && p == 1 
        seed = 1
    end
    if q == 7 
        seed = 5
    end
    metadata = joinpath(fpath,"princeton/data_w$(w0)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    push!(metadatas,metadata)
end
Δs= plot_spectra_collective(metadatas;savename="strain_s_1_t_3_spectrum.png")

fig = figure(figsize=(3.2,3))
plot(ϕs,Δs,"b-^")
xlim([0,0.55])
ylim([0,10])
ylabel("Δ (meV)")
xlabel(L"ϕ/ϕ_0")
tight_layout()
savefig("gap_vs_flux.pdf")
display(fig)
close(fig)
## BM basis 
seed = 1
p, q = 1 ,4
flag = "flavor"
νF = 2+ (2)*p/q
νstr = round(Int,1000*νF)
metadata = joinpath(fpath,"princeton/data_w07/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
plot_density_matrix_bm_valley_spin(metadata)
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
