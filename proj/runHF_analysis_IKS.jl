using PyPlot,JLD2
using Printf 
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF_IKS.jl"))
# include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))

#
# Info and folder name
# ------------------------------------------------------------------------------ # 
twist_angles = [105; collect(106:2:138)] 
twist_angle = 105
# for twist_angle in twist_angles
# dir = "/media/xiaoyuw@ad.magnet.fsu.edu/Data/Code/TBG_HartreeFock/"
# dir = "/Volumes/Data/Code/TBG_HartreeFock/zeeman/"
# dir1 = "/Volumes/Data/Reruns/"
# dir = "/Volumes/Data/Code/TBG_HartreeFock/MinHao/"
dir = ""
dir1 = ""
foldername = dir*"$(twist_angle)_strain"
# foldername = dir*"$(twist_angle)_strain"
params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)

# ----------------------------------Hartree Fock spectrum-------------------------------------------- # 
s,t = 0,-4
p,q = 1, 3
νF = (s)+(t)*p/q
νstr = round(Int,1000*νF)
metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="init_HF_$(p)_$(q)_nu_$(νstr)",_printinfo=true)

plot_spectra(metadata;savename="test.png",lines=[-10.0,4])
plot_spectrav3(metadata;savename="test.png",lines=[-2.5,2.5])
# plot_density_matrix_bm(metadata,ik=1)
# test_tL2_breaking(metadata)
# plot_density_matrix_global_order_parameters(metadata)

hf = load(metadata,"hf");
H = hf.H[(4q+1):6q,(4q+1):6q,:];
# H = hf.H[(1):2q,(1):2q,:];
U1 = zeros(ComplexF64,size(H,1),size(H,2),q*size(H,3))
tmpU1 = reshape(U1,size(H,1),size(H,2),q,size(H,3))
for ik in 1:size(H,3)
    F = eigen(Hermitian(H[:,:,ik]))
    for iq in 1:q
        tmpU1[:,:,iq,ik] = F.vectors
    end
end

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
mtg_data = "MinHao/NonInt/105_strain/_$(p)_$(q)_mtg_1_0_metadata.jld2";
function plot_realspace_cdw(metadata::String,mtg_data::String,ϵ0::Float64;γ::Float64=0.5)
    # γ = 1e3
    hf = load(metadata,"hf");
    U2 = zeros(ComplexF64,size(hf.H));
    energies = zeros(Float64,size(hf.ϵk))
    for ik in 1:size(hf.H,3) 
        F = eigen(Hermitian(view(hf.H,:,:,ik)))
        U2[:,:,ik] = F.vectors 
        energies[:,ik] = F.values
    end 
    ldim = size(load(dir1*hf.metadata[1],"Vec"),1)
    U1 = zeros(ComplexF64,ldim*hf.nη*hf.ns,size(hf.H,2),size(hf.H,3))
    U1tmp = reshape(U1,:,hf.nη,hf.ns,2hf.q,hf.nη,hf.ns,size(hf.H,3))
    for iη in 1:hf.nη 
        jldopen(dir1*hf.metadata[iη]) do file 
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
rvec, ldos  = plot_realspace_cdw(metadata,mtg_data,μ+4);
# mtg = load(mtg_data,"MTG");


# -------------------------------------------------------------------------------------------------------------------------- #

fig, ax = subplots(1,1,figsize=(5,8))
ldos1 = reshape(sum(ldos,dims=3),size(rvec))
# ldos1 = ldos[:,:,1]
# ldos2 = ldos[:,:,2]
ldos1 ./= maximum(ldos1)

for jj in -3:3
    pl=ax.pcolormesh(real(rvec).+jj*real(params.a1), imag(rvec).+jj*imag(params.a1),ldos1,vmin=0,vmax=1,cmap="bwr")
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



using DelimitedFiles

writedlm("Xgrid.txt",real(rvec))
writedlm("Ygrid.txt",imag(rvec))
writedlm("NormLDOS_1_6_-0.667_-3_-4meV_ConductionBand.txt",ldos1)
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
        F = eigen(Hermitian(view(hf.H0,:,:,ik)))
        U2[:,:,ik] = F.vectors 
        energies[:,ik] = F.values
    end 
    ldim = size(load(dir*hf.metadata[1],"Vec"),1)
    U1 = zeros(ComplexF64,ldim*hf.nη*hf.ns,size(hf.H,2),size(hf.H,3))
    U1tmp = reshape(U1,:,hf.nη,hf.ns,2hf.q,hf.nη,hf.ns,size(hf.H,3))
    for iη in 1:hf.nη 
        jldopen(dir*hf.metadata[iη]) do file 
            tmp = reshape(file["Vec"],ldim,2hf.q,hf.q,hf.nq^2)
            for is in 1:hf.ns 
                U1tmp[:,iη,is,:,iη,is,:] = view(tmp,:,:,1,:)
            end
        end
    end

    Utot = zeros(ComplexF64,size(U1))
    for ik in 1:size(U1,3)
        Utot[:,:,ik] = view(U1,:,:,ik) * view(U2,:,:,ik)
    end

    # Utot = zeros(ComplexF64,size(U2))
    # for ik in 1:size(U2,3)
    #     Utot[:,:,ik] = view(U2,:,:,ik)
    # end

    ldosLL = zeros(Float64,size(Utot,1))
    for ib in eachindex(ldosLL)
        ldosLL[ib] = sum(  abs2.(Utot[ib,:,:]) .* (sign.(energies .- ϵ1) .+1) .* (sign.(ϵ2 .- energies) .+1) ) / 4
    end
    return ldosLL ./ sum(ldosLL)
end

μ = load(metadata,"hf").μ
ldosLL  = plot_orbital_decomposition(metadata,-10.0,4.0);
ldosLL  = plot_orbital_decomposition(metadata,-2.5,0.0);


# ----------------------------------------------------------------
ldosLL_toplot = sum(reshape(ldosLL,:,p,2,2,2),dims=(2,3,4,5))[:];
# ldosLL_toplot = sum(reshape(ldosLL,:,2,2),dims=(2,3))[:];
nLL = 25*q÷p
idx = (-(nLL-1)):(nLL-1)
# idx = (-q+1):q
tmp = ldosLL_toplot[1:2:end]
tmp = [ldosLL_toplot[(end-1):-2:2]; tmp ]
writedlm("_$(p)_$(q)_orbital_decomposition_NonInt.txt",[idx tmp])


fig, ax = subplots(figsize=(4,3))
data1 = readdlm("_$(p)_$(q)_orbital_decomposition_NonInt.txt")
ax.plot(data1[:,1],data1[:,2]*100,"k-o",label="NonInt",ms=3)
data2 = readdlm("_$(p)_$(q)_orbital_decomposition_Int.txt")
ax.plot(data2[:,1],data2[:,2]*100,"b-o",label="(0,-4)",ms=3)
# ax.plot(data2[:,1],(data2[:,2]-data1[:,2])*100,"b-",label="diff",ms=2)
ax.legend()
ax.set_xlim([-15,15])
ax.set_xlabel("LL index")
ax.set_ylabel("Orbital weight (%)")
tight_layout()
savefig("_$(p)_$(q)_orbital_decomposition.png",dpi=600)
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


writedlm("-2_-2_HF_energy_diff_meta_stable_stable.txt",[data2 data1 ])

fig = figure(figsize=(4,3));
plot(data2,data1,"b-o",label="(-2,-1)"); 
legend()
xlabel(L"ϕ/ϕ_0")
ylabel(L"\rm ΔE_{HF}\ (meV)")
tight_layout()
display(fig); 
close(fig) 