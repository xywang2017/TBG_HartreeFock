using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
## Hartree Fock related 
params = Params()
initParamsWithStrain(params)

##
# w0s = ["00","02","03","05","06","07"]
w0s = ["00"]
w0snum = [0.0,0.2,0.3,0.5,0.6,0.7]
σz = []
p,q = 1, 8
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

fig = figure(figsize=(4,3)) 
plot(w0snum,reshape(σz,6,2)[:,1],"r-x",label="Chern")
plot(w0snum,reshape(σz,6,2)[:,2],"b-+",label="Flavor")
xlabel(L"w_0/w_1")
ylabel("⟨σzτz⟩")
title(L"ν=2.25,\ ϕ/ϕ_0=1/8")
legend()
tight_layout()
savefig("sublattice_polarization_v2.pdf",transparent=true)
display(fig)
close(fig)

fig = figure(figsize=(4,3))
data = zeros(Float64,length(w0s),3)
for ii in eachindex(w0s)
    w0 = w0s[ii]
    metadata1 = joinpath(fpath,"data_w$(w0)/_1_8/_chern_init_HF_1_8_nu_2250.jld2")
    metadata2 = joinpath(fpath,"data_w$(w0)/_1_8/_flavor_init_HF_1_8_nu_2250.jld2")
    data[ii,:] = [w0snum[ii] load(metadata1,"iter_energy")[end] load(metadata2,"iter_energy")[end]]
end
plot(data[:,1],data[:,2],"r-x",label="Chern")
plot(data[:,1],data[:,3],"b-+",label="Flavor")
xlabel(L"w_0/w_1")
ylabel("HF energy (meV)")
title(L"ν=2.25,\ ϕ/ϕ_0=1/8")
legend()
tight_layout()
savefig("hf_ground_state_quarter_flux_3250.pdf",transparent=true)
display(fig)
close(fig)


# HF energies at 1/4 and nu=1
data = [0.0 -51.29769545436033 -47.90256717821864;
        0.2 -50.890211134013704 -48.07554847054353;
        0.3 -50.46208073294238 -48.335005607404696;
        0.5 -49.57966993203583 -49.436512815452865;
        0.6 -49.15527663242192 -50.33236991559181;
        0.7 -48.239931777016544 -51.305786576313224]
# HF energies at 1/4 and nu=1.25
data = [0.0 -46.64867793444494 -48.79561815154276 ;
        0.2 -46.13247354504587 -48.636675207555434 ;
        0.3 -45.547385019860066 -48.4770812911432 ;
        0.5 -44.00941593088412 -48.22116249219488 ;
        0.6 -42.98546975215713 -48.187657184746435;
        0.7 -42.07981988464848 -48.099659491709446 ]
# HF energies at 1/4 and nu=1.75
data = [0.0 -44.51305735522851 -40.56474516094672;
        0.2 -43.73712628991553 -39.992833634780744;
        0.3 -42.8096500489483 -39.25946260858327;
        0.5 -40.07497806934608 -36.746957417805035;
        0.6 -38.21512222626289 -34.83998243479271;
        0.7 -35.59959989046577 -35.79403984619825]
# HF energies at 1/4 and nu=2.5
data = [0.0 -33.5070669178403 -30.30515943761121;
        0.2 -31.92974771397471 -29.153221740587295;
        0.3 -29.919898347502496 -27.62526564213171;
        0.5 -23.14169577753535 -22.031841740679308;
        0.6 -18.05606975306055 -17.529339501272585;
        0.7 -11.372540967050064 -11.5331885551484]

# HF energies at 1/4 and nu=3.25
data = [0.0 -18.27949695933558 -16.307055904871056;
        0.2 -15.454995404294394 -13.61605231734179;
        0.3 -11.759974096159215 -10.018986886743605;
        0.5 1.336259234069885 1.3367878377847142;
        0.6 11.492480191062779 11.492494531027477;
        0.7 24.640413834352643 24.641659227960055]

##
#plot everything
νs = collect(0.0:0.125:4.0)
νstrs = round.(Int,1000 .* νs)
fig = figure(figsize=(4,3)) 
for iν in eachindex(νs)
    ν = νs[iν]
    νstr = νstrs[iν] 
    metadata = joinpath(fpath,"data/_chern_init_HF_1_8_nu_$(νstr).jld2")
    if isfile(metadata)
        ϵk = load(metadata,"spectrum")
        chern = load(metadata,"chern")
        iter_err = load(metadata,"iter_err")[end]
        if iter_err < 1e-3
            scatter(ones(length(ϵk))*ν,sort(ϵk[:]),c=chern,cmap="Spectral_r",s=6,vmin=-1,vmax=1)
        end
    end
end
xlabel("ν")
ylabel("E (meV)")
# ylim([-55,70])
tight_layout()
# savefig("chern_spectrum.pdf",transparent=true)
display(fig)
close(fig)

# 
νs = collect(0.0:0.125:4.0)
νstrs = round.(Int,1000 .* νs)
hf_energies_chern = Float64[] 
hf_energies_flavor = Float64[]
fig = figure(figsize=(4,3)) 
for iν in eachindex(νs)
    ν = νs[iν]
    νstr = νstrs[iν] 
    metadata = joinpath(fpath,"data_chiral/_chern_init_HF_1_4_nu_$(νstr).jld2")
    push!(hf_energies_chern, load(metadata,"iter_energy")[end] )
    metadata = joinpath(fpath,"data_chiral/_flavor_init_HF_1_4_nu_$(νstr).jld2")
    push!(hf_energies_flavor, load(metadata,"iter_energy")[end] )
end
plot(νs,hf_energies_chern,"r-o",label="HF Chern",ms=4)
plot(νs,hf_energies_flavor,"b-o",label="HF Flavor",ms=4)
legend()
xlabel("ν")
ylabel("E (meV)")
tight_layout()
# savefig("flavor_chern_energy_comparison.pdf",transparent=true)
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

