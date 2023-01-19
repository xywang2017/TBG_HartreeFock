using PyPlot
using Printf
using DelimitedFiles
using JLD
fpath = joinpath(pwd(),"B0")
include(joinpath(fpath,"libs/HF_mod.jl"))

##
lk = 18
params = Params(ϵ=0.00,Da=0,dθ=1.05π/180,w1=96.056,w0=96.056*0.7,vf=2135.4)
initParamsWithStrain(params)
latt = Lattice()
initLattice(latt;lk=lk)
blk = HBM()
initHBM(blk,latt,params;lg=9)

kvec = reshape((real(latt.kvec)*params.g1 + imag(latt.kvec)*params.g2) ./ abs(params.g1),lk,lk)
# save(joinpath(fpath,"data/bareBM.jld"),"H0",blk.hbm)
ϵ0 = reshape(blk.hbm,blk.nη,blk.nb,lk,lk)

fig = figure(figsize=(4,3))
pl=contourf(real(kvec),imag(kvec),ϵ0[1,2,:,:],cmap="Spectral_r",levels=20)
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
for i in 1:2,j in 1:2
    color = (i==1) ? "r" : "b"
    plot(real(kvec[:,iΓ]),ϵ0[i,j,:,iΓ],".-",c=color)
    # plot(real(kcut),[ϵ[i,j,m,end-m+1] for m in 1:lk],".-",c=color)
end
ylabel(L"E")
# savefig("preHF.pdf",transparent=true)
display(fig)
close(fig)


### Hartree Fock part 
νs = collect(range(0,4,21))
hf = HartreeFock()
ν = 0.0
# for ν in νs
#     if (ν==1.2)
        run_HartreeFock(hf,blk,ν)
        save(joinpath(fpath,"data/nu$(ν).jld"),"H",hf.H,"P",hf.P,"ϵk",hf.ϵk)
#     end
# end

ν = 0.0
fp = joinpath(fpath,"data/nu$(ν).jld")
hf.ϵk = load(fp,"ϵk")
hf.H = load(fp,"H")
hf.P = load(fp,"P")
kvec = reshape((real(latt.kvec)*params.g1 + imag(latt.kvec)*params.g2) ./ abs(params.g1),lk,lk)
ϵ = reshape(hf.ϵk,8,lk,lk)

energies = range(minimum(hf.ϵk),maximum(hf.ϵk),2000)
Eν = 0.0
for iϵ in 1:length(energies)
    νrunning = sum((sign.(energies[iϵ] .- hf.ϵk).+1)./2) / blk.latt.nk .-4 
    if νrunning >= ν
        Eν = (energies[iϵ]+energies[iϵ-1])/2
        break
    end
end

Eν0 = 0.0
for iϵ in 1:length(energies)
    νrunning = 2*sum((sign.(energies[iϵ] .- blk.hbm).+1)./2) / blk.latt.nk .-4 
    if νrunning >= ν 
        Eν0 = (energies[iϵ]+energies[iϵ-1])/2
        break
    end
end

fig = figure(figsize=(4,3))
pcolormesh(real(kvec),imag(kvec),real(P_mix),cmap="coolwarm")
colorbar(label=L"⟨τ_xn_0s_0⟩")
cls = ["r","g","b","m"]
g1, g2, g3 = params.g1/abs(params.g1), params.g2/abs(params.g2), (params.g1+params.g2)/abs(params.g1)
lss = "solid"
for j in 5:1:8 
    contour(real(kvec),imag(kvec),reshape(hf.ϵk,blk.ns*blk.nfl,lk,lk)[j,:,:],levels=[Eν],colors=cls[j-4],linestyles=lss)
    # contour(real(kvec.+g1),imag(kvec.+g1),reshape(hf.ϵk,blk.ns*blk.nfl,lk,lk)[j,:,:],levels=[Eν],colors=cls[j-4],linestyles=lss)
    # contour(real(kvec.+g2),imag(kvec.+g2),reshape(hf.ϵk,blk.ns*blk.nfl,lk,lk)[j,:,:],levels=[Eν],colors=cls[j-4],linestyles=lss)
    # contour(real(kvec.+g3),imag(kvec.+g3),reshape(hf.ϵk,blk.ns*blk.nfl,lk,lk)[j,:,:],levels=[Eν],colors=cls[j-4],linestyles="solid")
end
lss = "dotted"
for i in 1:2, j in 1:2 
    if i==1
        contour(real(kvec),imag(kvec),ϵ0[i,j,:,:],levels=[Eν0],colors="m",linestyles=lss)
        # contour(real(kvec.+g1),imag(kvec.+g1),ϵ0[i,j,:,:],levels=[Eν0],colors="m",linestyles=lss)
        # contour(real(kvec.+g2),imag(kvec.+g2),ϵ0[i,j,:,:],levels=[Eν0],colors="m",linestyles=lss)
        # contour(real(kvec.+g3),imag(kvec.+g3),ϵ0[i,j,:,:],levels=[Eν0],colors="m",linestyles=lss)
    else
        contour(real(kvec),imag(kvec),ϵ0[i,j,:,:],levels=[Eν0],colors="g",linestyles=lss)
        # contour(real(kvec.+g1),imag(kvec.+g1),ϵ0[i,j,:,:],levels=[Eν0],colors="g",linestyles=lss)
        # contour(real(kvec.+g2),imag(kvec.+g2),ϵ0[i,j,:,:],levels=[Eν0],colors="g",linestyles=lss)
        # contour(real(kvec.+g3),imag(kvec.+g3),ϵ0[i,j,:,:],levels=[Eν0],colors="g",linestyles=lss)
    end
end
axis("equal")
display(fig)
close(fig)




###
s0 = ComplexF64[1 0; 0 1]
s1 = ComplexF64[0 1; 1 0]
s2 = ComplexF64[0 -1im; 1im 0]
s3 = ComplexF64[1 0; 0 -1]

P = reshape(hf.P,2,2,2,2,2,2,lk^2); 
P_I = reshape(diagm(ones(8)),2,2,2,2,2,2,1)

P_valley_s0 = 0.5 * reshape(sum((P.+0.5*P_I).*reshape(transpose(s0),1,2,1,1,2,1,1), dims=(2,5) )
            ,4,4,lk^2);
P_valley_s1 = 0.5 * reshape(sum(P.*reshape(transpose(s1),1,2,1,1,2,1,1), dims=(2,5) )
            ,4,4,lk^2);
P_valley_s2 = 0.5 * reshape(sum(P.*reshape(transpose(-1im*s2),1,2,1,1,2,1,1), dims=(2,5) )
            ,4,4,lk^2);
P_valley_s3 = 0.5 * reshape(sum((P.+0.5*P_I).*reshape(transpose(s3),1,2,1,1,2,1,1), dims=(2,5) )
            ,4,4,lk^2);

P_spin_s0 = 0.5 * reshape(sum((P.+0.5*P_I).*reshape(transpose(s0),2,1,1,2,1,1,1), dims=(1,4) )
            ,4,4,lk^2);
P_spin_s1 = 0.5 * reshape(sum(P.*reshape(transpose(s1),2,1,1,2,1,1,1), dims=(1,4) )
            ,4,4,lk^2);
P_spin_s2 = 0.5 * reshape(sum(P.*reshape(transpose(-1im*s2),2,1,1,2,1,1,1), dims=(1,4) )
            ,4,4,lk^2);
P_spin_s3 = 0.5 * reshape(sum(P.*reshape(transpose(s3),2,1,1,2,1,1,1), dims=(1,4) )
            ,4,4,lk^2);

P_band_s0 = 0.5 * reshape(sum(P.*reshape(transpose(s0),1,1,2,1,1,2,1), dims=(3,6) )
            ,4,4,lk^2);
P_band_s1 = 0.5 * reshape(sum(P.*reshape(transpose(s1),1,1,2,1,1,2,1), dims=(3,6) )
            ,4,4,lk^2);
P_band_s2 = 0.5 * reshape(sum(P.*reshape(transpose(-1im*s2),1,1,2,1,1,2,1), dims=(3,6) )
            ,4,4,lk^2);
P_band_s3 = 0.5 * reshape(sum(P.*reshape(transpose(s3),1,1,2,1,1,2,1), dims=(3,6) )
            ,4,4,lk^2);


P_valley = (0.5)^3* reshape(sum(reshape((P.+0.5P_I),8,8,:).*reshape(kron(s0,kron(s1,s0)),8,8,1),dims=(1,2) ),lk,lk)
P_band = (0.5)^3* reshape(sum(reshape((P.+0.5P_I),8,8,:).*reshape(kron(s0,kron(s0,s0)),8,8,1),dims=(1,2) ),lk,lk)
P_spin = (0.5)^3* reshape(sum(reshape((P.+0.5P_I),8,8,:).*reshape(kron(s0,kron(s0,s1)),8,8,1),dims=(1,2) ),lk,lk)

P_avg = reshape(sum(reshape(P,8,8,:),dims=3),8,8) / lk^2


P_mix = (0.5)^3* reshape(sum(reshape((P.+0.5P_I),8,8,:).*reshape(kron(s0,kron(s1,s0)),8,8,1),dims=(1,2) ),lk,lk);
println(norm(P_mix))

#######################################################
fig = figure(figsize=(4,3))
pl=pcolormesh(real(kvec),imag(kvec),ϵ[6,:,:],cmap="Spectral_r")
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

maximum(ϵ[4,:,:])
minimum(ϵ[5,:,:])
minimum(ϵ0[1,2,:,:])
maximum(ϵ0[1,1,:,:])

fig = figure(figsize=(4,3))
iΓ = (lk%2==0) ? (lk÷2) : ((lk-1)÷2+1)
xx = [collect(1:iΓ);iΓ.+collect(1:(iΓ-1))*sqrt(3)] # M -> Γ -> K -> M
# for i in 1:2,j in 1:2
#     ϵ0cut = [ϵ0[i,j,1:iΓ,iΓ]; [ϵ0[i,j,iΓ-m,iΓ+m] for m in 1:(iΓ-1)]]
#     plot(xx,ϵ0cut,"k.")
# end

for i in 1:8
    ϵcut = [ϵ[i,1:iΓ,iΓ]; [ϵ[i,iΓ-m,iΓ+m] for m in 1:(iΓ-1)]] ./V0
    plot(xx,ϵcut,"b.")
end
xticks([0.5,iΓ,iΓ+(iΓ-1+0.5)*sqrt(3)*2/3,iΓ+(iΓ-1+0.5)*sqrt(3)],["M","Γ","K","M"])
axvline(0.5,ls=":",c="gray")
axhline(0,ls=":",c="gray")
axvline(iΓ,ls=":",c="gray")
axvline(iΓ+(iΓ-1+0.5)*sqrt(3)*2/3,ls=":",c="gray")
axvline(iΓ+(iΓ-1+0.5)*sqrt(3),ls=":",c="gray")
ylabel(L"E")
# ylim([-80,80])
axhline(0.8)
tight_layout()
savefig("zerofield_HF_nu0.pdf",transparent=true)
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