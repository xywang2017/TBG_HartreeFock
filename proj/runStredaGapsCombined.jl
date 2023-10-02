using PyPlot,JLD2,Printf
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))

dir = "/media/xiaoyuw@ad.magnet.fsu.edu/Data/Code/TBG_HartreeFock/"

w0 = "07"
ϕs = sort(unique([p//q for q in 1:12 for p in 1:q]))
ϕs = ϕs[ϕs.<=0.5]
ϕs = ϕs[ϕs.>0.1]
sts = []
for s in -3:3, t in -12:12
    push!(sts,[s,t])
    push!(sts,[-s,-t])
end
sts = unique(sts)

# -------------------------Streda Line Plot ---------------------------------- # 
fig, ax = subplots(2,3,figsize=(10,6),sharex=true,sharey=true)
twist_angles = [138,132,128,124,120,105]
for i in eachindex(twist_angles)
    twist_angle = twist_angles[i]
    r, c = (i-1)÷3 + 1, (i-1)%3 + 1
    θ0 = @sprintf "%.2f" twist_angle*0.01
    # ax[r,c].set_title(L"θ=%$(θ0)^\circ")
    ax[r,c].text(-3.95,0.55,L"θ=%$(θ0)^\circ",size=11,c="k")
    foldername = dir*"zeeman/$(twist_angle)_nostrain"
    # for lines in -4:0
    #     ax[r,c].axvline(lines,ls=":",c="gray",lw=0.5)
    # end
    for ϕ in ϕs 
        p,q = numerator(ϕ), denominator(ϕ)
        gaps = Float64[]
        fillings = sort([st[1]+st[2]*p/q for st in sts])
        fillings = unique(round.(fillings,digits=10))
        fillings = fillings[fillings .<1e-5]
        ns = Float64[]
        for ν in fillings 
            νstr = round(Int,1000*ν)
            metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="nu_$(νstr)")
            if !isempty(metadata)
                push!(ns,ν)
                Δ,_=computegap(metadata;savename="test.png")
                push!(gaps,Δ)
                # if occursin("_tL_",metadata)
                #     println(metadata)
                # end
            end
        end
        ax[r,c].scatter(ns,ones(length(ns))*ϕ,s=gaps.^2/10,c="k",edgecolor="none")
    end
    if r==1 && c==1 
        ax[r,c].scatter([-1.0],[0.56],s=10^2/10,c="tab:red",edgecolor="none")
        ax[r,c].text(-0.9,0.55,"10 meV",size=11,c="tab:red")
    end
end

ax[1,1].set_xlim([-4.1,0.1])
ax[1,1].set_ylim([0.05,0.59])

for c in 1:3 
    ax[2,c].set_xlabel(L"n/n_s")
end
for r in 1:2
    ax[r,1].set_ylabel(L"ϕ/ϕ_0")
end
tight_layout()
savefig(joinpath(fpath,"fig1_nostrain.png"),transparent=false,dpi=600)
display(fig)
close(fig)


# ------------------------------- sublattice polarization analysis versus angle for a given flux --------------------------- 
twist_angles = [105; collect(106:2:138)]
# twist_angles = [105;120;124;128;132;138]
ϕ, s, t =1//8, -3, -1

# Pzs_strain = Float64[]
Pzs_bounds = ComplexF64[]
Pzs_symmetric = Float64[]
Pzs = Float64[]
Pzs_diag = Float64[]
Pzs_offdiag = Float64[]

for i in eachindex(twist_angles)
    twist_angle = twist_angles[i]
    p,q = numerator(ϕ), denominator(ϕ)
    foldername = dir*"zeeman/$(twist_angle)_strain"
    # ν = 0 + (-4)*ϕ
    # νstr = round(Int,1000*ν)
    # # metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="bm_cascade_tL_init_HF_$(p)_$(q)_nu_$(νstr)")
    # metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="$(p)_$(q)_nu_$(νstr)")
    # hf = load(metadata,"hf");
    # H = reshape(hf.H,8q,8q,:);
    # Σz = reshape(hf.Σz0,8q,8q,:);
    # nmax = round(Int,(ν+4)/8 *size(H,1)*size(H,3))
    # σz = 0.0 
    # for ik in 1:size(H,3)
    #     for i in 1:4
    #         vals, vec = eigen(Hermitian(H[((i-1)*(2q)+1):(i*(2q)),((i-1)*(2q)+1):(i*(2q)),ik]))
    #         σzτz = vec' * Σz[((i-1)*(2q)+1):(i*(2q)),((i-1)*(2q)+1):(i*(2q)),ik] * vec 
    #         σz  += sum(diag(σzτz)[1:(q-p)])
    #     end
    # end
    # push!(Pzs_symmetric,real(σz)/((q-p)*size(H,3)*4))

    
    ν = s + t*ϕ
    νstr = round(Int,1000*ν)
    metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="nu_$(νstr)")
    hf = load(metadata,"hf");
    P = reshape(hf.P,8q,8q,:);
    Σz = reshape(hf.Σz0,8q,8q,:);
    nmax = round(Int,(ν+4)/8 *size(P,1)*size(P,3))

    σz_upper = 0.0 
    for ik in 1:size(hf.P,3)
        tmpΣz =  eigvals(Hermitian(Σz[1:(2q),1:(2q),ik]))
        σz_upper += sum(tmpΣz[(q+p+1):(2q)])
    end

    σz_lower = 0.0 
    for ik in 1:size(hf.P,3)
        tmpΣz = diag(Σz[1:(2q),1:(2q),ik])
        σz_lower += sum(tmpΣz[1:(q-p)])
    end
    push!(Pzs_bounds,(σz_upper+1im*σz_lower)/((q-p)*size(P,3)))

    # σz = 0.0 
    # for ik in 1:size(hf.P,3)
    #     tmpP =  P[:,:,ik]+0.5I 
    #     tmpP = Diagonal(tmpP)
    #     σz += tr(transpose(tmpP) * view(Σz,:,:,ik)) / nmax 
    # end
    # push!(Pzs_diag,real(σz))

    # σz = 0.0 
    # for ik in 1:size(hf.P,3)
    #     tmpP =  P[:,:,ik]+0.5I 
    #     tmpP .= tmpP .- Diagonal(tmpP)
    #     σz += tr(transpose(tmpP) * view(Σz,:,:,ik)) / nmax 
    # end
    # push!(Pzs_offdiag,real(σz))

    σzmid = sum([tr(( transpose(view(P,:,:,ik))+0.5I ) * view(Σz,:,:,ik)) /nmax for ik in 1:size(P,3) ] ) 
    push!(Pzs,real(σzmid)) 
end


fig = figure(figsize=(3,2.5))
# fig = figure(figsize=(10,6))
plot(twist_angles.*0.01,Pzs,"bo",ms=3,label="($(s),$(t))")
plot(twist_angles.*0.01,real(Pzs_bounds),":",label="full Chern pol.")
plot(twist_angles.*0.01,imag(Pzs_bounds),":",label="BM valence")
# plot(twist_angles.*0.01,Pzs_symmetric,":",label="(0,-4) symm.")

legend(loc="upper right",fontsize=8)
xlabel("θ")
ylabel(L"\rm ⟨σ_zτ_z⟩")
tight_layout()
savefig("nostrain_Pz_s_$(s)_t_$(t).png",dpi=600,transparent=true)
display(fig)
close(fig)

