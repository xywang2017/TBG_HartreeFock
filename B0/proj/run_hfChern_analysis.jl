using PyPlot
using Printf
using JLD2
fpath = pwd() 
include(joinpath(fpath,"B0/libs/HFChern_mod.jl"))
include(joinpath(fpath,"B0/libs/plot_helpers.jl"))

twist_angle = 1.05
_is_strain = "strain"
foldername = @sprintf "%d_%s" round(Int,twist_angle*100) _is_strain
ν = 0.0
νstr = round(Int,1000*ν)
# ------------------ Specification ------------------ #
lk = 21
params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)
latt = Lattice()
initLattice(latt,params;lk=lk)


bm_path = joinpath(fpath,"$(foldername)/B0/bm_lk$(lk).jld2")
hf_path = find_lowest_energy_datafile("$(foldername)/B0";test_str="hf_$(νstr)_lk21",_printinfo=false)

# ----------------- Hartree-Fock dispersion part ---------------- # 
hf = load(hf_path,"hf");
println("HF energy: ", load(hf_path,"iter_energy")[end])
println("HF convergence: ",load(hf_path,"iter_err")[end])
kvec = reshape(latt.kvec ./ abs(params.g2),lk,lk)
# kvec = reshape(hf.latt.k1,:,1) .+ 1im*reshape(hf.latt.k2,1,:) 
ϵ0 = reshape(hf.ϵk,hf.nt,lk,lk)
# plot_contour_maps(kvec,ϵ0[4,:,:],points=ComplexF64[0+0im],contourlines=[1.],limits=Float64[])
# iΓ = (lk%2==0) ? (lk÷2) : ((lk-1)÷2+1)
# kcut = real(kvec[:,iΓ])
# Ecut = [ϵ0[j,i,iΓ] for j in 1:size(ϵ0,1),i in 1:size(ϵ0,3)];
# plot_energy_cuts(kcut,Ecut,lines=[hf.μ])

fig, ax = subplots(figsize=(6,4),subplot_kw=Dict("projection"=>"3d"))
ax.view_init(elev=20, azim=-50, roll=0)
pl = 0
for i in 1:2:8
    pl = ax.plot_surface(real(kvec),imag(kvec),ϵ0[i,:,:],cmap="coolwarm",vmin=minimum(ϵ0),vmax=maximum(ϵ0))
end
# ax.set_xticks(collect(-0.4:0.2:0.4))
# ax.set_yticks(collect(-0.4:0.2:0.4))
# ax.set_xlim([-0.55,0.55])
# ax.set_ylim([-0.55,0.55])
ax.set_xlabel(L"k_x / L_2")
ax.set_ylabel(L"k_y /L_2")
ax.set_zlabel("E (meV)")
colorbar(pl,shrink=0.4,location="left")
tight_layout()
savefig("tmp.png",dpi=500,transparent=true)
display(fig)
close(fig)



### all the chemical potentials 
νs = collect(-4.0:0.125:4.0)
μs = Float64[]
actual_νs = Float64[]
for ν in νs 
    νstr = round(Int,1000*ν)
    hf_path = find_lowest_energy_datafile("$(foldername)/B0";test_str="hf_$(νstr)_",_printinfo=false)
    if ispath(hf_path)
        hf = load(hf_path,"hf");
        push!(actual_νs,round(Int,(hf.ν+4)/8*size(hf.H,1)*size(hf.H,3))/(size(hf.H,1)*size(hf.H,3))*8 - 4)
        push!(μs,hf.μ)
    end
end
fig = figure(figsize=(4,3))
plot(sort(actual_νs),μs[sortperm(actual_νs)],"b^",ms=3,lw=1)
xlabel(L"n/n_s")
ylabel("μ (meV)")
# axhline(0,ls=":",c="gray")
# axvline(0,ls=":",c="gray")
xticks(collect(-3:3))
ylim([-30,30])
yticks(collect(-20:10:20))
# ylim([0,14])
xlim([-4.1,4.1])
tight_layout()
savefig("cascade_strain_gapless.pdf",dpi=500,transparent=true)
display(fig)
close(fig)

