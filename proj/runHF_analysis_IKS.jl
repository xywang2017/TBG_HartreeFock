using PyPlot,JLD2
using Printf 
fpath = pwd()
# include(joinpath(fpath,"libs/MagneticFieldHF_IKS.jl"))
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))

#
# Info and folder name
# ------------------------------------------------------------------------------ # 
twist_angles = [105; collect(106:2:138)] 
twist_angle = 120
# for twist_angle in twist_angles
dir = "/media/xiaoyuw@ad.magnet.fsu.edu/Data/Code/TBG_HartreeFock/zeeman/"
dir = "/Volumes/Data/Code/TBG_HartreeFock/zeeman/"
dir = ""
foldername = dir*"$(twist_angle)_nostrain"
params = Params(ϵ=0.00,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)

# ----------------------------------Hartree Fock spectrum-------------------------------------------- # 
s,t = -2,0
p,q = 3, 8
νF = (s)+(t)*p/q
νstr = round(Int,1000*νF)
metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="init_HF_$(p)_$(q)_nu_$(νstr)",_printinfo=true)

# plot_spectra(metadata;savename="test.png")
plot_density_matrix_bm(metadata,ik=1)
# test_tL2_breaking(metadata)
plot_density_matrix_global_order_parameters(metadata)

# ----------------------------------IKS Analysis-------------------------------------------- # 
nq = 12÷q 
if q == 7 
    nq = 2  
end
q1s = collect(0:(nq*q-1))
energies = []
for q1 in q1s 
    metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="_$(q1)_0_random_init_HF_$(p)_$(q)_nu_$(νstr)",_printinfo=true)
    push!(energies,load(metadata,"iter_energy")[end])
end
fig = figure(figsize=(4,3))
plot(q1s./(nq*q),energies .- energies[1],"k--^")
xlabel(L"q_1")
ylabel(L"E(q_1)-E(0)\ (meV)")
xlim([-0.09,1.09])
tight_layout()
# savefig("test.png",dpi=600)
display(fig)
close(fig)


P = reshape(load(metadata,"hf").P,2q,4,2q,4,:);

p_subblock = P[:,3,:,2,:];

_pratio = p_subblock[:,:,2]./ p_subblock[:,:,1]

angle.(_pratio) /(2π) /(1/8)

# ------------------------------ IKS energy vs wavevector 
qs = reshape(0:7,:,1) .+ 1im*reshape(0:7,1,:)
energies = zeros(Float64,size(qs))
for iq in eachindex(qs)
    q1, q2 = real(qs[iq]), imag(qs[iq])
    metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="_$(q1)_0_$(q2)_flavor_tL_init_HF_$(p)_$(q)_nu_$(νstr)",_printinfo=true)
    energies[iq] = load(metadata,"iter_energy")[end]
end

fig = figure(figsize=(3.6,3))
pl = pcolor(real(qs), imag(qs),energies,cmap="bwr")
colorbar(pl,shrink=0.8)
xlabel(L"q_1")
ylabel(L"q_2")
tight_layout()
display(fig)
close(fig)

# ---------------------------- real space density modulation at a given energy ---------------------- # 
mtg_data = "NonInt/120_strain/_$(p)_$(q)_mtg_0_0_metadata.jld2";
function plot_realspace_cdw(metadata::String,mtg_data::String,ϵ0::Float64;γ::Float64=0.5)
    hf = load(metadata,"hf");
    U2 = zeros(ComplexF64,size(hf.H));
    energies = zeros(Float64,size(hf.ϵk))
    for ik in 1:size(hf.H,3) 
        F = eigen(Hermitian(view(hf.H,:,:,ik)))
        U2[:,:,ik] = F.vectors 
        energies[:,ik] = F.values
    end 
    ldim = size(load(dir*hf.metadata[1],"Vec"),1)
    U1 = zeros(ComplexF64,ldim*hf.nη*hf.ns,size(hf.H,2),size(hf.H,3))
    U1tmp = reshape(U1,:,hf.nη,hf.ns,2hf.q,hf.nη,hf.ns,size(hf.H,3))
    for iη in 1:hf.nη 
        jldopen(dir*hf.metadata[iη]) do file 
            for is in 1:hf.ns 
                U1tmp[:,iη,is,:,iη,is,:] = reshape(file["Vec"],ldim,2hf.q,:)
            end
        end
    end

    Utot = zeros(ComplexF64,size(U1))
    for ik in 1:size(U1,3)
        Utot[:,:,ik] = view(U1,:,:,ik) * view(U2,:,:,ik)
    end

    # real space info 
    mtg = load(mtg_data,"MTG");
    uvec = reshape(Utot,:,size(hf.H,2),size(hf.H,3));
    ψ = zeros(ComplexF64,2,mtg.coord.nr,size(hf.H,2),size(hf.H,3));  #layerxsublatticexvalleyxspin , total bands, momentum, position
    W = reshape(mtg.W,2,:,size(hf.H,3),mtg.coord.nr,hf.nη,1);

    for ik in 1:size(ψ,4), ib in 1:size(ψ,3), ir in 1:mtg.coord.nr, sublatt in 1:2 
        ψ[sublatt,ir,ib,ik] = sum(view(W,sublatt,:,ik,ir,:,:).*reshape(view(uvec,:,ib,ik),:,hf.nη,hf.ns))
    end

    rvec = reshape(mtg.coord.z,mtg.coord.lr,:)
    # rvec = reshape(mtg.coord.x1,:,1) .+ 1im*reshape(mtg.coord.x2,1,:)
    ir = reshape(1:length(rvec),size(rvec))
    ldos = zeros(Float64,size(rvec,1),size(rvec,2),2)
    for x in 1:size(ldos,1), y in 1:size(ldos,2), jj in 1:2 
        ldos[x,y,jj] = sum(  abs2.(ψ[jj,ir[x,y],:,:]) ./ ((ϵ0 .- energies).^2 .+ γ^2) )* γ/π
        # ldos[x,y,jj] = sum(  abs2.(ψ[jj,ir[x,y],:,:]) .* (sign.(hf.μ .- energies) .+1 ) ) / 2
    end
    return rvec, ldos
end

μ = load(metadata,"hf").μ
rvec, ldos  = plot_realspace_cdw(metadata,mtg_data,0.0);
# mtg = load(mtg_data,"MTG");


# -------------------------------------------------------------------------------------------------------------------------- #

fig, ax = subplots(1,1,figsize=(5,8))
ldos1 = reshape(sum(ldos,dims=3),size(rvec))
# ldos1 = ldos[:,:,2]
# ldos2 = ldos[:,:,2]
ldos1 ./= maximum(ldos1)

for jj in -3:3
    pl=ax.pcolormesh(real(rvec).+jj*real(params.a1), imag(rvec).+jj*imag(params.a1),ldos1, cmap="bwr",vmin=0,vmax=1)
end
points = reshape(-3:3,:,1)*params.a1 .+ reshape(-3:3,1,:)*params.a2 
ax.scatter(real(points),imag(points),s=5,c="k")
ax.axis("equal")
ax.axis("off")
# colorbar(pl,shrink=0.6)
tight_layout()
savefig("test.png",dpi=600,transparent=true)
display(fig)
close(fig)





# ------------
hf = load(metadata,"hf");
fig = figure(figsize=(4,3))
kvec = reshape(hf.latt.k1[1:8],:,1) * hf.params.g1 .+ reshape(hf.latt.k2[1:8],1,:)*hf.params.g2
contourf(real(kvec),imag(kvec),reshape(hf.ϵk[9,:],8,8))
# pcolormesh(real(kvec),imag(kvec),reshape(hf.ϵk[8,:],8,8))
colorbar()
# contour(real(kvec),imag(kvec),reshape(hf.ϵk[8,:],8,8),[hf.μ])

axis("equal")
display(fig)
savefig("abc.png",dpi=500)
close(fig)



# ----------------------------------------------------LL orbital decomposition for states in a given energy interval ------------------------------- #
function plot_orbital_decomposition(metadata::String,ϵ1::Float64,ϵ2::Float64;γ::Float64=0.5)
    hf = load(metadata,"hf");
    U2 = zeros(ComplexF64,size(hf.H));
    energies = zeros(Float64,size(hf.ϵk))
    for ik in 1:size(hf.H,3) 
        F = eigen(Hermitian(view(hf.H,:,:,ik)))
        U2[:,:,ik] = F.vectors 
        energies[:,ik] = F.values
    end 
    ldim = size(load(dir*hf.metadata[1],"Vec"),1)
    U1 = zeros(ComplexF64,ldim*hf.nη*hf.ns,size(hf.H,2),size(hf.H,3))
    U1tmp = reshape(U1,:,hf.nη,hf.ns,2hf.q,hf.nη,hf.ns,size(hf.H,3))
    for iη in 1:hf.nη 
        jldopen(dir*hf.metadata[iη]) do file 
            for is in 1:hf.ns 
                U1tmp[:,iη,is,:,iη,is,:] = reshape(file["Vec"],ldim,2hf.q,:)
            end
        end
    end

    # Utot = zeros(ComplexF64,size(U1))
    # for ik in 1:size(U1,3)
    #     Utot[:,:,ik] = view(U1,:,:,ik) * view(U2,:,:,ik)
    # end

    Utot = zeros(ComplexF64,size(U2))
    for ik in 1:size(U2,3)
        Utot[:,:,ik] = view(U2,:,:,ik)
    end

    ldosLL = zeros(Float64,size(Utot,1))
    for ib in eachindex(ldosLL)
        ldosLL[ib] = sum(  abs2.(Utot[ib,:,:]) .* (sign.(energies .- ϵ1) .+1) .* (sign.(ϵ2 .- energies) .+1) ) / 4
    end
    return ldosLL ./ sum(ldosLL)
end

μ = load(metadata,"hf").μ
ldosLL  = plot_orbital_decomposition(metadata,-8.0,0.0);


# ----------------------------------------------------------------
# ldosLL_toplot = sum(reshape(ldosLL,:,p,2,2,2),dims=(2,3,4,5))[:];
ldosLL_toplot = sum(reshape(ldosLL,:,2,2),dims=(2,3))[:];
fig, ax = subplots(figsize=(4,3))
# nLL = 25*q÷p
# idx = (-(nLL-1)):(nLL-1)
idx = (-q+1):q
# tmp = ldosLL_toplot[1:2:end] + ldosLL_toplot[1:2:end]
# tmp = [ldosLL_toplot[(end-1):-2:2] + ldosLL_toplot[(end-1):-2:2]; tmp ]
# ax.plot(idx,tmp*100,"g.-",ms=2)
ax.plot(idx,ldosLL_toplot,"g.-",ms=2)
# ax.set_xlim([-79,79])
# ax.set_ylim([-2,19])
# ax.set_yticks(collect(0:5:15))
ax.set_xlabel("LL index")
ax.set_ylabel("Orbital weight (%)")
ax.set_title("1.05 strain (0,-2) $(p)/$(q)")
tight_layout()
# savefig("105_flux_$(p)_$(q).png",dpi=600)
display(fig)
close(fig)


# ------------------------------------------------------------------ #
ldosLL_toplot = reshape(ldosLL,:,2,2,2);
fig, ax = subplots(2,2,figsize=(5,4),sharex=true,sharey=true)
lblstr = ["K↑" "K↓";"K'↑" "K'↓";]
nLL = 25*q÷p
for η in 1:2, s in 1:2 
    idx = (-(nLL-1)):(nLL-1)
    tmp = ldosLL_toplot[1:2:end,1,η,s] + ldosLL_toplot[1:2:end,2,η,s]
    tmp = [ldosLL_toplot[(end-1):-2:2,1,η,s] + ldosLL_toplot[(end-1):-2:2,2,η,s]; tmp ]
    # ax[η,s].plot(0:(nLL-1),ldosLL_toplot[1:2:end,2,η,s]+ldosLL_toplot[1:2:end,1,η,s],"g.",ms=2,label=lblstr[η,s])
    # ax[η,s].plot(-1:-1:-(nLL-1),ldosLL_toplot[2:2:end,2,η,s]+ldosLL_toplot[2:2:end,1,η,s],"g.",ms=2)
    ax[η,s].plot(idx,tmp*100,"g-",ms=2,label=lblstr[η,s])
    # ax[η,s].set_yscale("log")
    if η ==2 
        ax[η,s].set_xlabel("LL index")
    end
    if s == 1
        ax[η,s].set_ylabel("Orbital weight (%)")
    end
    ax[η,s].legend()
    ax[1,1].set_title("1.05 strain (0,-2) $(p)/$(q)")
    # ax[η,s].set_xlim(-60,60)
end
tight_layout()
savefig("105_flux_$(p)_$(q).png",dpi=600)
display(fig)
close(fig)
