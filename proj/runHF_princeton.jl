using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
## Hartree Fock related 
params = Params()
initParamsWithStrain(params)

##

w0cs = Float64[]
w0s = ["07"]
w0snum = [0.7]
σz = []
p,q = 1, 4
νF = 2 + 2*p/q
νstr = round(Int,1000*νF)
for w0 in w0s
    metadata = joinpath(fpath,"princeton/old/data_w$(w0)/_$(p)_$(q)/_flavor_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    println(load(metadata,"iter_energy")[end])
    println(load(metadata,"iter_err")[end])
    # plot_hf_iterations(metadata)
    ϵk = load(metadata,"spectrum")
    σzτz = load(metadata,"chern")
    fig = figure(figsize=(4,3)) 
    scatter(ones(length(ϵk)),ϵk[:],c=σzτz[:],s=6,cmap="coolwarm",vmin=-1,vmax=1)
    ϵsorted = sort(ϵk[:])
    ν = eachindex(ϵsorted) ./ length(ϵsorted)
    i = 1
    while (νF+4)/8 > ν[i]
        i += 1
    end
    ϵF = (ϵsorted[i+1] + ϵsorted[i])/2 
    Δ = (ϵsorted[i+1] - ϵsorted[i]) 
    println(i," ",length(ϵsorted))
    println(Δ)
    axhline(ϵF)
    display(fig)
    close(fig)
    # push!(σz, plot_spectra(ϵk,σzτz,νF,params;savename=joinpath(fpath,"princeton/figures/chern_spectrum_$(p)_$(q)_$(w0).png")))
end 

for w0 in w0s
    metadata = joinpath(fpath,"princeton/data_w$(w0)/_$(p)_$(q)/_flavor_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    println(load(metadata,"iter_energy")[end])
    println(load(metadata,"iter_err")[end])
    # plot_hf_iterations(metadata)
    ϵk = load(metadata,"spectrum")
    σzτz = load(metadata,"chern")
    push!(σz, plot_spectra(ϵk,σzτz,νF,params;savename=joinpath(fpath,"princeton/figures/flavor_spectrum_$(p)_$(q)_$(w0).png")))
end 

Σz = readdlm("princeton/spectrum_gap.txt")
ϕs = 1 ./ [collect(4:8);10]
fig = figure(figsize=(4,3))
pcolormesh(ϕs,w0snum,Σz[1:6,:])
colorbar()
xlim([0,0.3])
ylim([0,0.8])
ylabel(L"w_0/w_1")
xlabel(L"ϕ/ϕ_0")
tight_layout()
savefig("princeton/figures/spectrum_gap_chern.pdf",transparent=true)
display(fig)
close(fig)

fig = figure(figsize=(4,3))
for i in eachindex(w0s)
    str = w0s[i]
    if i  == 1
        plot(ϕs,Σz[i,:],"rp-",ms=4,label="c$(str)")
        plot(ϕs,Σz[i+6,:],"rp-",ms=4,label="f$(str)")
    elseif i==6
        plot(ϕs,Σz[i,:],"b^-",ms=4,label="c$(str)")
        plot(ϕs,Σz[i+6,:],"b^-",ms=4,label="f$(str)")
    else
        plot(ϕs,Σz[i,:],".--",ms=2,label="c$(str)")
        plot(ϕs,Σz[i+6,:],".--",ms=2,label="f$(str)")
    end
end
legend()
xlim([0,0.3])
ylabel(L"Δ (E_c)")
xlabel(L"ϕ/ϕ_0")
tight_layout()
savefig("princeton/figures/spectrum_gap_cuts.pdf",transparent=true)
display(fig)
close(fig)


fig = figure(figsize=(4,3)) 
plot(w0snum,reshape(σz,6,2)[:,1],"r-x",label="Chern")
plot(w0snum,reshape(σz,6,2)[:,2],"b-+",label="Flavor")
xlabel(L"w_0/w_1")
# ylabel("⟨σzτz⟩")
ylabel(L"Δ (E_c)")
title(L"s=2, t=2, ϕ/ϕ_0=%$(p)/%$(q)")
legend()
tight_layout()
# savefig("princeton/figures/sublattice_polarization_$(p)_$(q).pdf",transparent=true)
savefig("princeton/figures/spectrum_gap_$(p)_$(q).pdf",transparent=true)
display(fig)
close(fig)

### Hartree-Fock ground state energy 
fig = figure(figsize=(4,3))
data = zeros(Float64,length(w0s),3)
for ii in eachindex(w0s)
    w0 = w0s[ii]
    metadata1 = joinpath(fpath,"princeton/data_w$(w0)/_$(p)_$(q)/_chern_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    metadata2 = joinpath(fpath,"princeton/data_w$(w0)/_$(p)_$(q)/_flavor_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    data[ii,:] = [w0snum[ii] load(metadata1,"iter_energy")[end] load(metadata2,"iter_energy")[end]]
end
plot(data[:,1],data[:,2],"r-x",label="Chern")
plot(data[:,1],data[:,3],"b-+",label="Flavor")
xlabel(L"w_0/w_1")
ylabel("HF energy (meV)")
title(L"s=2, t=2, ϕ/ϕ_0=%$(p)/%$(q)")
legend()
tight_layout()
savefig("princeton/figures/hf_ground_state_$(p)_$(q).pdf",transparent=true)
display(fig)
close(fig)

## crossing point via linear interpolation
cnt = 1
for i in 1:size(data,1)
    if data[i,2] > data[i,3]
        cnt = i
        break 
    end
end
x1,y1 = data[cnt-1,1] ,data[cnt-1,2] - data[cnt-1,3]
x2,y2 = data[cnt,1], data[cnt,2] - data[cnt,3]
k = (y2-y1)/(x2-x1)
b = y1 
push!(w0cs, x1-b/k)

fig = figure(figsize=(4,3))
tmp = readdlm("princeton/w0cs.txt")
plot(tmp[:,1],tmp[:,2],"bx-")
ylabel(L"w_0/w_1")
xlabel(L"ϕ/ϕ_0")
ylim([0,0.8])
xlim([0,0.3])
tight_layout()
savefig("princeton/figures/first_order_phase_transition.pdf",transparent=true)
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
H0 = load(joinpath(fpath,"princeton/data_w06/_1_5/_flavor_init_HF_1_5_nu_0.jld2"),"H");
P = load(joinpath(fpath,"princeton/data_w06/_1_5/_flavor_init_HF_1_5_nu_400.jld2"),"P");
ik = 5
tmpH0 = reshape(view(H0,:,:,ik),10,4,10,4)
P0 = reshape(view(P,:,:,ik)+0.5I,10,4,10,4)
Pstrong = zeros(ComplexF64,10,10,4)
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

