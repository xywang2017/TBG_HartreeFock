using PyPlot
using JLD2
fpath = joinpath(pwd(),"B0")
include(joinpath(fpath,"libs/HF_mod.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))

prefix =1
# νs = collect(0.0:0.2:4.0)
ν = 2.0
νstr = round(Int,1000*ν)
# ------------------ Specification ------------------ #
lk = 19
params = Params(ϵ=0.002,Da=-4100,dθ=1.06π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)
latt = Lattice()
initLattice(latt,params;lk=lk)

bm_path = joinpath(fpath,"data/strain_bm_lk$(lk).jld2")
hf_path = joinpath(fpath,"data/$(prefix)_strain_hf_$(νstr)_lk$(lk).jld2")
# hf_path = "typical_starting_point.jld2"


# ----------------- valley-spin polarization info ----------------- # 
fig = figure(figsize=(2,10))
plot(hf.Δ,eachindex(hf.Δ),"b^")
# yticks(collect(eachindex(hf.Δ))[2:end],hf.Δstr[2:end])
yticks(collect(1:length(hf.Δ)))
axvline(0,c="gray")
xlim(-0.4,0.4)
tight_layout()
display(fig)
close(fig)

##### ------ 
ρ = hf.Δ[1]
α = hf.Δ[13]
Δ0 = hf.Δ[2:4]
Δ3 = hf.Δ[14:16]
ΔK = Δ0 + Δ3 
ΔKprime = Δ0 - Δ3
ρK = ρ+α
ρKprime = ρ-α