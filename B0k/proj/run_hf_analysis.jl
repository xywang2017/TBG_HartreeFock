using PyPlot
using JLD2
fpath = joinpath(pwd(),"B0k")
include(joinpath(fpath,"libs/HF_mod.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))

prefix =1
# νs = collect(0.0:0.2:4.0)
ν = 1.0
νstr = round(Int,1000*ν)
# ------------------ Specification ------------------ #
lk = 8
params = Params(ϵ=0.003,Da=-4100,dθ=1.06π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)
latt = Lattice()
initLattice(latt,params;lk=lk)

bm_path = joinpath(fpath,"data/strain_bm_lk$(lk).jld2")
hf_path = joinpath(fpath,"data/$(prefix)_strain_hf_$(νstr)_lk$(lk).jld2")
# hf_path = "typical_starting_point.jld2"

# ----------------- Hartree-Fock analysis part ---------------- # 
hf = load(hf_path,"hf");
plot_energies(hf.ϵk,hf.σzτz,lines=[hf.μ])


# ----------------- Hartree-Fock statistics ---------------- # 
fig = figure(figsize=(5,4)) 
pl = imshow(abs.(hf.P-Diagonal(hf.P)),cmap="Blues",origin="lower")
colorbar(pl)
axis("equal")
tight_layout()
display(fig)
close(fig)