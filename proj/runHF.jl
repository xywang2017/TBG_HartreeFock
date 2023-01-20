using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
## Hartree Fock related 
params = Params()
initParamsWithStrain(params)

##
w0s = ["00","02","03","05","06","07"]
w0snum = [0.0,0.2,0.3,0.5,0.6,0.7]
σz = []
p,q = 1, 4
νF = 2 + 2*p/q
νstr = round(Int,1000*νF)
for w0 in w0s
    metadata = joinpath(fpath,"data_w$(w0)/_$(p)_$(q)/_chern_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    println(load(metadata,"iter_energy")[end])
    println(load(metadata,"iter_err")[end])
    # plot_hf_iterations(metadata)
    ϵk = load(metadata,"spectrum")
    σzτz = load(metadata,"chern")
    push!(σz, plot_spectra(ϵk,σzτz,νF,params;savename=joinpath(fpath,"figures/chern_spectrum_$(w0).png")))
end 

for w0 in w0s
    metadata = joinpath(fpath,"data_w$(w0)/_$(p)_$(q)/_flavor_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    println(load(metadata,"iter_energy")[end])
    println(load(metadata,"iter_err")[end])
    # plot_hf_iterations(metadata)
    ϵk = load(metadata,"spectrum")
    σzτz = load(metadata,"chern")
    push!(σz, plot_spectra(ϵk,σzτz,νF,params;savename=joinpath(fpath,"figures/chern_spectrum_$(w0).png")))
end 

fig = figure(figsize=(4,3)) 
plot(w0snum,reshape(σz,6,2)[:,1],"r-x",label="Chern")
plot(w0snum,reshape(σz,6,2)[:,2],"b-+",label="Flavor")
xlabel(L"w_0/w_1")
ylabel("⟨σzτz⟩")
title(L"s=2, t=2, ϕ/ϕ_0=%$(p)/%$(q)")
legend()
tight_layout()
savefig("sublattice_polarization_$(p)_$(q).pdf",transparent=true)
display(fig)
close(fig)

fig = figure(figsize=(4,3))
data = zeros(Float64,length(w0s),3)
for ii in eachindex(w0s)
    w0 = w0s[ii]
    metadata1 = joinpath(fpath,"data_w$(w0)/_$(p)_$(q)/_chern_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    metadata2 = joinpath(fpath,"data_w$(w0)/_$(p)_$(q)/_flavor_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    data[ii,:] = [w0snum[ii] load(metadata1,"iter_energy")[end] load(metadata2,"iter_energy")[end]]
end
plot(data[:,1],data[:,2],"r-x",label="Chern")
plot(data[:,1],data[:,3],"b-+",label="Flavor")
xlabel(L"w_0/w_1")
ylabel("HF energy (meV)")
title(L"s=2, t=2, ϕ/ϕ_0=%$(p)/%$(q)")
legend()
tight_layout()
savefig("hf_ground_state_$(p)_$(q).pdf",transparent=true)
display(fig)
close(fig)

## plot density matrix 
P = load(metadata,"P");
ik = 1
fig = figure(figsize=(4,3))
pl=imshow(abs.(P[:,:,ik] + 0.5I),origin="lower")
colorbar(pl)
tight_layout()
display(fig)
close(fig)

## strong coupling basis 
H0=load(joinpath(fpath,"data_w07/_flavor_init_HF_1_4_nu_250.jld2"),"H");
ik = 8
tmpH0 = reshape(view(H0,:,:,ik),8,4,8,4)
P0 = reshape(view(P,:,:,ik)+0.5I,8,4,8,4)
Pstrong = zeros(ComplexF64,8,8,4)
for iηs in 1:4
    F = eigen(Hermitian(tmpH0[:,iηs,:,iηs]))
    vec = F.vectors
    Pstrong[:,:,iηs] = vec' * P0[:,iηs,:,iηs] * vec 
end 
fig,ax = subplots(2,2,figsize=(8,8))
states = ["K↑","K'↑","K↓","K'↓"]
for r in 1:2, c in 1:2 
    pl=ax[r,c].imshow(abs.(Pstrong[:,:,r+2(c-1)]),vmin=0,vmax=1)
    colorbar(pl,ax=ax[r,c],fraction=0.046, pad=0.04)
    ax[r,c].set_title(states[r+2(c-1)])
end
tight_layout()
display(fig)
close(fig)

