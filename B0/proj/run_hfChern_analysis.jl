using PyPlot
using JLD2
fpath = pwd() 
include(joinpath(fpath,"B0/libs/HFChern_mod.jl"))
include(joinpath(fpath,"B0/libs/plot_helpers.jl"))

prefix = 1
flag = "random"
phi = 0
strain = 2
# νs = collect(0.0:0.2:4.0)
ν = -2.0
νstr = round(Int,1000*ν)
# ------------------ Specification ------------------ #
lk = 16
# params = Params(ϵ=0.00,Da=0,dθ=1.06π/180,w1=110,w0=77,vf=2482)
params = Params(ϵ=0.001*strain,Da=-4100,dθ=1.05π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)
latt = Lattice()
initLattice(latt,params;lk=lk)

bm_path = joinpath(fpath,"feldman/B0/data/strain$(strain)/phi$(phi)/bm_lk$(lk).jld2")
hf_path = joinpath(fpath,"feldman/B0/data/strain$(strain)/phi$(phi)/$(prefix)_$(flag)_hf_$(νstr)_lk$(lk).jld2")

# ----------------- Hartree-Fock dispersion part ---------------- # 
hf = load(hf_path,"hf");
iter_energy = load(hf_path,"iter_energy");
println(iter_energy[end])
# kvec = reshape(latt.kvec ./ abs(params.g1),lk,lk)
kvec = reshape(hf.latt.k1,:,1) .+ 1im*reshape(hf.latt.k2,1,:) 
ϵ0 = reshape(hf.ϵk,hf.nt,lk,lk)
# plot_contour_maps(kvec,ϵ0[4,:,:],points=[params.Kt/abs(params.g1)],contourlines=[hf.μ],limits=Float64[])
iΓ = (lk%2==0) ? (lk÷2) : ((lk-1)÷2+1)
kcut = real(kvec[:,iΓ])
Ecut = ϵ0[:,:,iΓ]
plot_energy_cuts(kcut,Ecut,lines=[hf.μ])

# ----------------- valley-spin-chern polarization info ----------------- # 
s0 = ComplexF64[1 0;0 1]
s1 = ComplexF64[0 1;1 0]
s2 = ComplexF64[0 -1im;1im 0]
s3 = ComplexF64[1 0;0 -1]
paulis = [s1,s2,s3,s0]
# Δ = zeros(size(hf.ϵk))
δs = zeros(Float64,4,4,4,size(hf.ϵk,2))
for ik in 1:size(hf.ϵk,2)
    # F = eigen(Hermitian(view(hf.H,:,:,ik)))
    # Δ[:,ik] = real(diag(F.vectors'*kron(s,kron(s1,s0))*F.vectors))
    # Δ[:,ik] = sqrt.(real(diag(F.vectors'*kron(s0,kron(s1,s0))*F.vectors)).^2 + 
    #             real(diag(F.vectors'*kron(s0,kron(s2,s0))*F.vectors)).^2 )
    # check_Hermitian(hf.P[:,:,ik])
    for i in 1:4, j in 1:4, k in 1:4
        δs[i,j,k,ik] = real(tr(transpose(hf.P[:,:,ik]+I/2)*kron(paulis[i],kron(paulis[j],paulis[k]))))/8
    end
    # tmpP = zeros(ComplexF64,8,8)
    # # for i in 1:4, j in 1:4, k in 1:4
    # #     tmpP .+= δs[i,j,k,ik]* kron(paulis[i],kron(paulis[j],paulis[k]))
    # # end
    # for i in 1:4, j in 1:4
    #     if !(i==4 && j==4)
    #         tmpP .+= δs[i,j,4,ik]* kron(paulis[i],kron(paulis[j],paulis[4]))
    #     end
    # end
    # term1 = zeros(ComplexF64,8,8)
    # for i in 1:3 
    #     term1 .+= δs[4,4,i,ik]*kron(paulis[4],kron(paulis[4],paulis[i]))
    # end
    # tmpP .= (δs[4,4,4,ik]-1)*I + (I+tmpP) * (I+term1)
    # if norm(tmpP .- transpose(hf.P[:,:,ik]+I/2) ) > 1e-10
    #     println(norm(abs.(tmpP .- transpose(hf.P[:,:,ik]+I/2) )))
    # end
end
# Δ .= hf.σzτz
fig = figure(figsize=(3,3))
pl = pcolor(real(kvec),imag(kvec),reshape(δs[4,3,4,:].*δs[2,4,4,:],hf.latt.lk,hf.latt.lk),cmap="bwr")
colorbar(pl)
display(fig)
close(fig)
# plot_energy_cuts_with_order_parameters(kcut,Ecut,
#                 reshape(Δ,:,lk,lk)[:,:,iΓ],lines=[hf.μ])

# plot_contour_maps(kvec,reshape(Δ,:,lk,lk)[1,:,:],points=[params.Kt/abs(params.g1)],contourlines=[hf.μ])
# plot_density_maps(kvec,collect(reshape(δs[2,4,4,:],lk,lk)),
#             points=[params.Kt/abs(params.g1)],contourlines=Float64[],limits=Float64[])
plot_density_maps_collective(kvec,collect(reshape(δs,4,4,4,lk,lk)),
            points=[params.Kt/abs(params.g1)],contourlines=Float64[],limits=Float64[])

plot_density_maps_collectivev0(kvec,collect(reshape(δs,4,4,4,lk,lk)),
            points=[params.Kt/abs(params.g1)],contourlines=Float64[],limits=Float64[])

fig = figure(figsize=(4,3))
quiver(real(kvec),imag(kvec),collect(reshape(δs[1,4,4,:],lk,lk)),collect(reshape(δs[2,4,4,:],lk,lk)))
axis("equal")
tight_layout()
savefig("test.pdf")
display(fig)
close(fig)
# ----------------- order parameters ----------------- # 
fig, ax = subplots(1,2,sharex=true,figsize=(4,5))
s = length(hf.Δ)
ax[1].axvline(0,c="gray")
ax[1].plot(hf.Δ[2:s÷2],eachindex(hf.Δ)[2:s÷2],"b^")
ax[1].set_yticks(collect(eachindex(hf.Δ))[2:s÷2])
ax[1].set_yticklabels(hf.Δstr[2:s÷2])
ax[2].axvline(0,c="gray")
ax[2].plot(hf.Δ[(s÷2+1):end],eachindex(hf.Δ)[(s÷2+1):end],"b^")
ax[2].set_yticks(collect(eachindex(hf.Δ))[(s÷2+1):end])
ax[2].set_yticklabels(hf.Δstr[(s÷2+1):end])
# ax[2].set_xlim(-0.2,0.2)
tight_layout()
savefig("test.pdf")
display(fig)
close(fig)

# ---------------------------- Hartree only ---------------------------------- # 
Δ = zeros(size(hf.ϵk))
δ = zeros(size(hf.ϵk,2))
ϵk = zeros(size(hf.ϵk))
for ik in 1:size(hf.ϵk,2)
    tmpH = reshape(view(hf.H,:,:,ik),4,2,4,2)
    for j in 1:4
        ϵk[(2j-1):(2j),ik] = eigvals(Hermitian(tmpH[j,:,j,:])) 
    end
end
# Δ .= hf.σzτz

# plot_energy_cuts_with_order_parameters(kcut,Ecut,
#                 reshape(Δ,:,lk,lk)[:,:,iΓ],lines=[hf.μ])

plot_contour_maps(kvec,reshape(ϵk,:,lk,lk)[6,:,:],points=[params.Kt/abs(params.g1)],contourlines=[hf.μ])
# ----------------- Hartree-Fock statistics ---------------- # 
iter_oda = load(hf_path,"iter_oda");
iter_err = load(hf_path,"iter_err");
iter_energy = load(hf_path,"iter_energy");
println(iter_energy[end])
println(iter_err[end])

fig,ax = subplots(figsize=(4,3))
ax.plot(eachindex(iter_err),iter_err,"b-",label="err",ms=2)
ax1 = ax.twinx()
ax1.plot(eachindex(iter_energy),iter_energy,"r-",label="oda",ms=3)
ax.set_ylabel("Iter Err")
ax.set_yscale("log")
ax1.set_ylabel("Iter Energy")
# ax1.set_yscale("log")
tight_layout()
savefig("test.pdf")
display(fig)
close(fig)

### all the chemical potentials 
νs = collect(-3.9:0.1:3.9)
lk = 20
phi = 0
strain = 3
μs = Float64[]
actual_νs = Float64[]
for ν in νs 
    νstr = round(Int,1000*ν)
    hf_path = joinpath(fpath,"feldman/B0/data/strain$(strain)/phi$(phi)/1_random_hf_$(νstr)_lk$(lk).jld2")
    if ispath(hf_path)
        hf = load(hf_path,"hf");
        push!(actual_νs,round(Int,(hf.ν+4)/8*size(hf.H,1)*size(hf.H,3))/(size(hf.H,1)*size(hf.H,3))*8 - 4)
        push!(μs,hf.μ)
    end
end
fig = figure(figsize=(5,3))
for i in -4:4 
    axvline(i,ls=":",c="gray")
end
axhline(0,ls=":",c="gray")
plot(sort(actual_νs),μs[sortperm(actual_νs)],"b-o",ms=2)
xlabel("ν")
ylabel("μ")
axhline(0,ls="--",c="gray")
# ylim([0,30])
# ylim([0,14])
xlim([-4.2,4.2])
tight_layout()
savefig("cascade_strain_phi$(phi).pdf")
display(fig)
close(fig)