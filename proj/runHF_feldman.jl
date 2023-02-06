using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
## Hartree Fock related 
params = Params()
initParamsWithStrain(params)

##
flag = "flavor"
w0cs = Float64[]
w0s = ["06"]
w0snum = [0.6]
σz = []
p,q = 1, 4
νF = 0 + (2)*p/q
νstr = round(Int,1000*νF)
for w0 in w0s
    metadata = joinpath(fpath,"feldman/data_w$(w0)/_$(p)_$(q)/_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    println(load(metadata,"iter_energy")[end])
    println(load(metadata,"iter_err")[end])
    # plot_hf_iterations(metadata)
    ϵk = load(metadata,"spectrum")
    σzτz = load(metadata,"chern")
    push!(σz, plot_spectra(ϵk,σzτz,νF,params;savename=joinpath(fpath,"feldman/figures/$(flag)_spectrum_$(p)_$(q)_$(w0).png")))
end 

## BM basis 
P = load(joinpath(fpath,"feldman/data_w06/_$(p)_$(q)/_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2"),"P");
ik = 3
P0 = reshape(view(P,:,:,ik)+0.5I,2q,4,2q,4);
fig,ax = subplots(2,2,figsize=(8,8))
states = ["K↑","K'↑","K↓","K'↓"]
for r in 1:2, c in 1:2 
    pl=ax[r,c].imshow(abs.(P0[:,r+2(c-1),:,r+2(c-1)]),vmin=0,vmax=1)
    colorbar(pl,ax=ax[r,c],fraction=0.046, pad=0.04)
    ax[r,c].set_title(states[r+2(c-1)])
end
tight_layout()
display(fig)
close(fig)

fig = figure(figsize=(6,6))
pl = imshow(abs.(reshape(P0,8q,:)),vmin=0,vmax=1)
colorbar(pl)
tight_layout()
display(fig)
close(fig)

## strong coupling basis
H0 = load(joinpath(fpath,"feldman/data_w06/_$(p)_$(q)/_$(flag)_init_HF_$(p)_$(q)_nu_0.jld2"),"H");
P = load(joinpath(fpath,"feldman/data_w06/_$(p)_$(q)/_$(flag)_init_HF_$(p)_$(q)_nu_500.jld2"),"P");
ik = 1
tmpH0 = view(H0,:,:,ik)
P0 = view(P,:,:,ik)+0.5I
F = eigen(Hermitian(tmpH0))
Pstrong = F.vectors' * P0 * F.vectors 

fig = figure(figsize=(6,6))
pl=imshow(abs.(Pstrong),vmin=0,vmax=1)
colorbar(pl,fraction=0.046, pad=0.04)
tight_layout()
display(fig)
close(fig)



## strong coupling basis valley spin resolvecd
H0 = load(joinpath(fpath,"feldman/data_w06/_1_5/_$(flag)_init_HF_1_5_nu_0.jld2"),"H");
P = load(joinpath(fpath,"feldman/data_w06/_1_5/_$(flag)_init_HF_1_5_nu_400.jld2"),"P");
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
