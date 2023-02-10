using PyPlot
using JLD2
fpath = joinpath(pwd(),"B0")
include(joinpath(fpath,"libs/HF_mod.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))

prefix = 2
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

# ----------------- Hartree-Fock analysis part ---------------- # 
hf = load(hf_path,"hf");
kvec = reshape(latt.kvec ./ abs(params.g1),lk,lk)
ϵ0 = reshape(hf.ϵk,hf.nt,lk,lk)
plot_contour_maps(kvec,ϵ0[3,:,:],points=[params.Kt/abs(params.g1)])
iΓ = (lk%2==0) ? (lk÷2) : ((lk-1)÷2+1)
kcut = real(kvec[:,iΓ])
Ecut = ϵ0[:,:,iΓ]
μ = find_chemicalpotential(hf.ϵk[:],(hf.ν+4)/8)
plot_energy_cuts(kcut,Ecut,lines=[μ])


# ----------------- Hartree-Fock statistics ---------------- # 
iter_oda = load(hf_path,"iter_oda");
iter_err = load(hf_path,"iter_err");
iter_energy = load(hf_path,"iter_energy");

fig,ax = subplots(figsize=(4,3))
ax.plot(eachindex(iter_err),iter_err,"b.",label="err",ms=2)
ax1 = ax.twinx()
ax1.plot(eachindex(iter_energy),iter_energy,"r-",label="oda",ms=3)
ax.set_ylabel("Iter Err")
ax.set_yscale("log")
ax1.set_ylabel("Iter Energy")
tight_layout()
savefig("test.pdf")
display(fig)
close(fig)