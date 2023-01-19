using PyPlot
using Printf
using DelimitedFiles
using JLD
using Random
fpath = joinpath(pwd(),"B0")
include(joinpath(fpath,"libs/HFvalleyK_mod.jl"))

##
lk = 20
params = Params(ϵ=0.002,Da=-4100,dθ=1.38π/180,w1=110.0,w0=88.0,vf=2680.0)
initParamsWithStrain(params)
latt = Lattice()
initLattice(latt;lk=lk)
blk = HBM()
initHBM(blk,latt,params;lg=9)

kvec = reshape((real(latt.kvec)*params.g1 + imag(latt.kvec)*params.g2) ./ abs(params.g1),lk,lk)
# save(joinpath(fpath,"data/bareBM.jld"),"H0",blk.hbm)
ϵ0 = reshape(blk.hbm,blk.nb,lk,lk)

fig = figure(figsize=(4,3))
pl=contourf(real(kvec),imag(kvec),ϵ0[2,:,:],cmap="Spectral_r",levels=20)
plot([real(params.Kt)]/abs(params.g1),[imag(params.Kt)]/abs(params.g1),"k+")
colorbar(pl)
xlabel(L"k_x")
ylabel(L"k_y")
axis("equal")
display(fig)
close(fig)


fig = figure(figsize=(4,3))
kcut = [kvec[i,end-i+1] for i in 1:lk]
iΓ = (lk%2==0) ? (lk÷2) : ((lk-1)÷2+1)
for i in 1:2
    color = (i==1) ? "r" : "b"
    plot(real(kvec[:,iΓ]),ϵ0[i,:,iΓ],".-",c=color)
    # plot(real(kcut),[ϵ[i,j,m,end-m+1] for m in 1:lk],".-",c=color)
end
ylabel(L"E")
# savefig("preHF.pdf",transparent=true)
display(fig)
close(fig)


### Hartree Fock part 
νs = collect(range(0,4,21))
hf = HartreeFock()
ν = 0.13
# for ν in νs
#     if (ν==1.2)
        run_HartreeFock(hf,blk,ν)
        save(joinpath(fpath,"datavalleyK/nu$(ν).jld"),"H",hf.H,"P",hf.P,"ϵk",hf.ϵk)
#     end
# end

ν = 0.13
fp = joinpath(fpath,"datavalleyK/nu$(ν).jld")
hf.ϵk = load(fp,"ϵk")
hf.H = load(fp,"H")
hf.P = load(fp,"P")
kvec = reshape((real(latt.kvec)*params.g1 + imag(latt.kvec)*params.g2) ./ abs(params.g1),lk,lk)
ϵ = reshape(hf.ϵk,2,lk,lk)
energies = range(minimum(hf.ϵk),maximum(hf.ϵk),2000)
Eν = 0.0
for iϵ in 1:length(energies)
    νrunning = sum((sign.(energies[iϵ] .- hf.ϵk).+1)./2) / blk.latt.nk .-1 
    if νrunning >= ν
        Eν = (energies[iϵ]+energies[iϵ-1])/2
        break
    end
end

Eν0 = 0.0
for iϵ in 1:length(energies)
    νrunning = sum((sign.(energies[iϵ] .- blk.hbm).+1)./2) / blk.latt.nk .-1
    if νrunning >= ν
        Eν0 = (energies[iϵ]+energies[iϵ-1])/2
        break
    end
end

fig = figure(figsize=(4,3))
cls = ["r","b","g","m"]
g1, g2, g3 = params.g1/abs(params.g1), params.g2/abs(params.g2), (params.g1+params.g2)/abs(params.g1)
lss = "solid"
for j in 1:2 
    contour(real(kvec),imag(kvec),reshape(hf.ϵk,2,lk,lk)[j,:,:],levels=[Eν],colors=cls[j],linestyles=lss)
    contour(real(kvec.+g1),imag(kvec.+g1),reshape(hf.ϵk,2,lk,lk)[j,:,:],levels=[Eν],colors=cls[j],linestyles=lss)
    contour(real(kvec.+g2),imag(kvec.+g2),reshape(hf.ϵk,2,lk,lk)[j,:,:],levels=[Eν],colors=cls[j],linestyles=lss)
    contour(real(kvec.+g3),imag(kvec.+g3),reshape(hf.ϵk,2,lk,lk)[j,:,:],levels=[Eν],colors=cls[j],linestyles="solid")
end
lss = "solid"
for j in 1:2 
    contour(real(kvec),imag(kvec),ϵ0[j,:,:],levels=[Eν0],colors="k",linestyles=lss)
    contour(real(kvec.+g1),imag(kvec.+g1),ϵ0[j,:,:],levels=[Eν0],colors="k",linestyles=lss)
    contour(real(kvec.+g2),imag(kvec.+g2),ϵ0[j,:,:],levels=[Eν0],colors="k",linestyles=lss)
    contour(real(kvec.+g3),imag(kvec.+g3),ϵ0[j,:,:],levels=[Eν0],colors="k",linestyles=lss)
end
axis("equal")
display(fig)
close(fig)

#######################################################
fig = figure(figsize=(4,3))
pl=pcolormesh(real(kvec),imag(kvec),ϵ[1,:,:],cmap="Spectral_r")
plot([real(params.Kt)]/abs(params.g1),[imag(params.Kt)]/abs(params.g1),"k+")
colorbar(pl)
xlabel(L"k_x")
ylabel(L"k_y")
axis("equal")
tight_layout()
display(fig)
close(fig)

Lm = sqrt(abs(params.a1)*abs(params.a2))*2.46e-10 
ee = 1.6e-19
ϵϵ = 8.8541878128e−12	
ϵr = 5.0
V0 = ee^2/(4π*ϵϵ*ϵr* Lm) /ee *1e3

ee^2/(4π*ϵϵ* abs(params.a1'*params.a2)*(2.46e-10)) /(1.6e-19) *1e3

fig = figure(figsize=(4,3))
iΓ = (lk%2==0) ? (lk÷2) : ((lk-1)÷2+1)
xx = [collect(1:iΓ);iΓ.+collect(1:(iΓ-1))*sqrt(3)] # M -> Γ -> K -> M
for j in 1:2
    ϵ0cut = [ϵ0[j,1:iΓ,iΓ]; [ϵ0[j,iΓ-m,iΓ+m] for m in 1:(iΓ-1)]]
    plot(xx,ϵ0cut,"k--")
end

for i in 1:2
    ϵcut = [ϵ[i,1:iΓ,iΓ]; [ϵ[i,iΓ-m,iΓ+m] for m in 1:(iΓ-1)]]
    plot(xx,ϵcut,"b-")
end
xticks([0.5,iΓ,iΓ+(iΓ-1+0.5)*sqrt(3)*2/3,iΓ+(iΓ-1+0.5)*sqrt(3)],["M","Γ","K","M"])
axvline(0.5,ls=":",c="gray")
axvline(iΓ,ls=":",c="gray")
axvline(iΓ+(iΓ-1+0.5)*sqrt(3)*2/3,ls=":",c="gray")
axvline(iΓ+(iΓ-1+0.5)*sqrt(3),ls=":",c="gray")
ylabel(L"E")
# ylim([-80,80])
tight_layout()
# savefig("unstrained_HF_nu0.pdf",transparent=true)
display(fig)
close(fig)


### H analysis 
tmpH = reshape(hf.H,2,2,2,2,2,2,blk.latt.nk);

# test valley coherence 
tmpH_ivc = norm(tmpH[:,1,:,:,2,:,:]) +  norm(tmpH[:,2,:,:,1,:,:])
ϵsK = zeros(Float64,4,blk.latt.nk)
for ik in 1:blk.latt.nk 
    ϵsK[:,ik] = real(eigvals(reshape(tmpH[:,1,:,:,1,:,ik],4,4))) 
end
ϵsK = reshape(ϵsK,4,lk,lk)

fig = figure(figsize=(4,3))
pl=contourf(real(kvec),imag(kvec),ϵsK[3,:,:],cmap="Spectral_r",levels=20)
plot([real(params.Kt)]/abs(params.g1),[imag(params.Kt)]/abs(params.g1),"k+")
colorbar(pl)
xlabel(L"k_x")
ylabel(L"k_y")
axis("equal")
display(fig)
close(fig)


fig = figure(figsize=(4,3))
iΓ = (lk%2==0) ? (lk÷2) : ((lk-1)÷2+1)
xx = [collect(1:iΓ);iΓ.+collect(1:(iΓ-1))*sqrt(3)] # M -> Γ -> K -> M
for i in 1:2,j in 1:2
    ϵ0cut = [ϵ0[i,j,1:iΓ,iΓ]; [ϵ0[i,j,iΓ-m,iΓ+m] for m in 1:(iΓ-1)]]
    plot(xx,ϵ0cut,"k--")
end

for i in 1:4
    ϵcut = [ϵsK[i,1:iΓ,iΓ]; [ϵsK[i,iΓ-m,iΓ+m] for m in 1:(iΓ-1)]]
    plot(xx,ϵcut,"b-")
end
xticks([0.5,iΓ,iΓ+(iΓ-1+0.5)*sqrt(3)*2/3,iΓ+(iΓ-1+0.5)*sqrt(3)],["M","Γ","K","M"])
axvline(0.5,ls=":",c="gray")
axvline(iΓ,ls=":",c="gray")
axvline(iΓ+(iΓ-1+0.5)*sqrt(3)*2/3,ls=":",c="gray")
axvline(iΓ+(iΓ-1+0.5)*sqrt(3),ls=":",c="gray")
ylabel(L"E")
# ylim([-80,80])
tight_layout()
# savefig("unstrained_HF_nu0.pdf",transparent=true)
display(fig)
close(fig)


### dos 
energies = collect(range(0,80,200))
γ = 2.0
dos0 = zeros(Float64,length(energies))
dos = zeros(Float64,length(energies))
ν0s = zeros(Float64,length(energies))
νs = zeros(Float64,length(energies))
for ii in eachindex(energies)
    dos0[ii] = 2*γ/π * sum(1 ./(γ^2 .+ (energies[ii] .- ϵ0[1,2,:,:]).^2))
    dos[ii] = γ/π * sum(1 ./(γ^2 .+ (energies[ii] .- ϵsK[3:4,:,:]).^2))
    ν0s[ii] = 2*sum((sign.(energies[ii] .- ϵ0[1,2,:,:]).+1)./2) / latt.nk
    νs[ii] = sum((sign.(energies[ii] .- ϵsK[3:4,:,:]).+1)./2) / latt.nk
end

fig = figure() 
plot(energies,dos0,"k--",label="non-int")
plot(energies,dos,"b-",label="hf")
xlabel("E (meV)")
ylabel("DoS")
legend()
tight_layout()
display(fig)
close(fig)


fig = figure() 
plot(2 .*ν0s,dos0,"k--",label="non-int")
plot(2 .*νs,dos,"b-",label="hf")
xlabel("ν")
ylabel("DoS")
legend()
tight_layout()
display(fig)
close(fig)