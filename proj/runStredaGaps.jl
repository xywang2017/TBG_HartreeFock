using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
# dir = "/media/xiaoyuw@ad.magnet.fsu.edu/Data/Code/TBG_HartreeFock/"
dir = "/Volumes/Data/Code/TBG_HartreeFock/zeeman/"
# dir = ""
# Info and folder name
# ------------------------------------------------------------------------------ # 
twist_angle = 120
foldername = dir*"$(twist_angle)_nostrain"
fname1 = dir*"MinHao/$(twist_angle)_strain"
# params = Params(ϵ=0.002,Da=0.0,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2125.6)
params = Params(ϵ=0.000,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)


w0 = "07"

ϕs = sort(unique([p//q for q in 1:12 for p in 1:q]))
ϕs = ϕs[ϕs.<=0.5]
# ϕs = ϕs[ϕs.>=1//8]
sts = []
for s in -3:3, t in -12:12
    push!(sts,[s,t])
end
sts = unique(sts)
μs = ComplexF64[]
# -------------------------Streda Line Plot ---------------------------------- # 
cs = ["r";"g";"b";"c";"m";"darkviolet";"tab:blue";
        "magenta";"peru";"tab:purple";"tab:olive";"deepskyblue";"seagreen";"gray"]
fig, ax = subplots(figsize=(5,4))
# for lines in -4:4
#     axvline(lines,ls=":",c="gray",lw=0.5)
# end
for iϕ in eachindex(ϕs)
    ϕ = ϕs[iϕ] 
    p,q = numerator(ϕ), denominator(ϕ)
    gaps = Float64[]
    fillings = sort([st[1]+st[2]*p/q for st in sts])
    fillings = unique(round.(fillings,digits=8))
    fillings = fillings[fillings .<1e-5]
    fillings = fillings[fillings .>-4]
    # fillings = fillings[abs.(fillings) .<4]
    fillings0 = round.(fillings,digits=3)
    ns = Float64[]
    for ν in fillings
        νstr = round(Int,1000*ν)
        if q <=12
            metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="nu_$(νstr).jld2")
        else
            metadata = find_lowest_energy_datafile("$(fname1)/_$(p)_$(q)";test_str="nu_$(νstr).jld2")
        end
        # println(metadata)
        if !isempty(metadata)
            push!(ns,ν)
            Δ,Pz=computegap(metadata;savename="test.png")
            push!(gaps,Δ)
            if ϕ == 1//12 
                push!(μs,ν+1im*load(metadata,"hf").μ)
            end
        end
        
    end
    # gaps[gaps .< 5.5] .= 0.0 
    ax.scatter(ns,ones(length(ns))*ϕ,s=gaps.^2 ./10,c="tab:blue",edgecolor="none")
    # if ! isempty(ns)
    #     writedlm("105_fillings_gaps/$(p)_$(q)_fillings_gaps.txt",[ns gaps])
    # end
end
ax.set_xlim([-4.3,0.3])
ax.set_ylim([-0.01,0.55])
ax.set_xlabel(L"n/n_s",fontsize=13)
ax.set_ylabel(L"ϕ/ϕ_0",fontsize=13)
ax2 = ax.twinx()
mn, mx = ax.get_ylim()
ax2.set_ylim(flux_conversion(mn,params), flux_conversion(mx,params))
ax2.set_yticks(collect(5:5:flux_conversion(mx,params)))
ax2.set_ylabel("B (T)",fontsize=13)
# scatter([-0.5-3*1/3],[1/3])
tight_layout()
savefig("$(twist_angle).png",transparent=false,dpi=600)
# savefig(joinpath(fpath,"$(foldername)/streda_line.png"),transparent=false,dpi=600)
display(fig)
close(fig)

# -----

fig = figure(figsize=(4,3))
plot(real(μs),imag(μs),"b-o",ms=3)
xlabel(L"n/n_s")
ylabel("μ (meV)")
tight_layout()
display(fig)
close(fig)
# -----------------------------Hofstadter spectrum plot ---------------------------- # 
Δss = []
energies1 = zeros(Float64,length(ϕs),3)
energies = zeros(Float64,length(ϕs),3)
sts = [[-1,-3],[-2,-2],[-3,-1]]
sts = [[-3,-1]]
for i in eachindex(sts) 
    st = sts[i]
    s,t = st[1], st[2]
    metadatas = String[]
    for iϕ in eachindex(ϕs)
        ϕ = ϕs[iϕ] 
        p,q = numerator(ϕ), denominator(ϕ)
        νstr = round(Int,1000*(s+t*p/q))
        if abs(s+t*p/q) < 4
            metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="tL_init_HF_$(p)_$(q)_nu_$(νstr)")
            energies1[iϕ,i] = load(metadata,"iter_energy")[end]
            metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="init_HF_$(p)_$(q)_nu_$(νstr)")
            energies[iϕ,i] = load(metadata,"iter_energy")[end]
            if !isempty(metadata)
                push!(metadatas,metadata)
                # println(load(metadata,"iter_energy")[end])
            end
        end
    end
    idx = Int[]
    # idx = collect(1:14)
    # idx = collect(1:length(ϕs))
    Δs= plot_spectra_collectivev2(metadatas;savename="spectrum_s$(s)_t$(t).png",titlestr="(s,t)=($(s),$(t))",indices=idx);
    push!(Δss,Δs)
end


# fig = figure(figsize=(2.5,2.5))
# colors = ["tab:blue","tab:red","tab:green","tab:orange","cyan"]
# for i in eachindex(Δss)
#     s, t = sts[i][1], sts[i][2]
#     plot(ϕs,Δss[i],"o",c=colors[i],ms=4,markeredgecolor="none",label="($(s),$(t))")
# end
# legend(fontsize=8,loc="upper left")
# xlabel(L"ϕ/ϕ_0")
# ylabel("Δ (meV)")
# xlim([-0.1,0.55])
# ylim([0,24])
# tight_layout()
# savefig("tmp_$(twist_angle).png",transparent=true,dpi=600)
# display(fig)
# close(fig)

fig,ax = subplots(2,2,figsize=(3,2.5),sharex=true,sharey=true)
colors = ["tab:blue","tab:red","tab:green","tab:orange","cyan"]
colors = ["k","k","k","k","k"]
for i in eachindex(Δss)
    s, t = sts[i][1], sts[i][2]
    ax[i].plot(ϕs,Δss[i],"o",c=colors[i],ms=3,markeredgecolor="none",label="($(s),$(t))")
    ax[i].legend(fontsize=8,loc="upper left")
    if (i-1)%2 +1 == 2 
        ax[i].set_xlabel(L"\rm ϕ/ϕ_0")
    end
    if (i-1)÷2 +1 == 1
        ax[i].set_ylabel("Δ (meV)")
    end
end
subplots_adjust(hspace=0,wspace=0)
ax[1].set_xlim([0,0.55])
ax[1].set_xticks([0.2,0.4])
ax[1].set_ylim([0,29])
# tight_layout()
savefig("tmp_$(twist_angle).png",transparent=true,dpi=600,bbox_inches="tight")
display(fig)
close(fig)


# ---------------- non interacting hofstadter spectrum weighted with a given Streda line density matrix
# plot spectrum 
function plot_LL_spectrum(params::Params;angle::Int=120)
    foldername0 =dir*"NonInt/$(angle)_strain"
    fig,ax = subplots(figsize=(2.5,2.5))
    # fig,ax = subplots(figsize=(6,4))
    # ϕs = [1//5]
    strs = ["K","Kprime"]
    pl = 0
    for ϕ in ϕs
        p, q = numerator(ϕ), denominator(ϕ)
        strs = ["K","Kprime"]
        colors = ["b","r"]
        tmp = []
        for iη in 1:1
            str = strs[iη]
            fname0 = joinpath(fpath,"$(foldername0)/_$(p)_$(q)_$(str)_metadata.jld2")
            energies = load(fname0,"E");
            Σz = load(fname0,"PΣz");
            weights = zeros(Float64,2q,size(energies,2),size(energies,3))
            for i2 in 1:size(energies,3), i1 in 1:size(energies,2)
                weights[:,i1,i2] = real(diag(Σz[:,:,i1,i2]))
            end
            energies = reshape(energies,2q,:)
            pl = ax.scatter(ones(length(energies[1:(q-p),:]))*ϕ,energies[1:(q-p),:],marker="o",s=1.5,edgecolor="none",c=[[0.8,0,0]])
            pl = ax.scatter(ones(length(energies[(q-p+1):end,:]))*ϕ,energies[(q-p+1):end,:],marker="o",s=1.5,edgecolor="none",c="gray")
            # pl = ax.scatter(ones(length(energies[:]))*ϕ,energies[:],marker="o",s=6,edgecolor="none",c=weights[:]*(3-2iη),vmin=-1,vmax=1,cmap="coolwarm")
        end
        
    end
    # colorbar(pl,fraction=0.06, pad=0.05)
    ax.set_ylabel("E (meV)")
    ax.set_xlabel(L"ϕ/ϕ_0")
    ax.set_xlim([0,0.55])
    ax.set_ylim([-45,45])
    tight_layout()
    savefig("tmp_$(angle).png",dpi=600,transparent=true)
    display(fig)
    close(fig)
    return nothing
end

plot_LL_spectrum(params,angle=132)



# ------------------- 
data = readdlm("Bc_theta.txt")
yvals = zeros(Float64,size(data,1),size(data,2)-1)
for i in 1:size(data,1)
    twist_angle = data[i,1] * 100
    params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
    initParamsWithStrain(params)
    for j in 2:size(data,2) 
        yvals[i,j-1] = flux_conversion(data[i,j],params)
    end
end

fig, ax = subplots(1,3,figsize=(4,2.5),sharex=true,sharey=true)
lblstr = ["(-3,-1)","(-2,-2)","(-1,-3)"]
colors = ["tab:blue","tab:red","tab:green","tab:orange","cyan"]
for j in 1:3
    ax[j].plot(data[:,1],yvals[:,j],"--x",label=lblstr[j],c=colors[j+1])
    ax[j].legend(loc="upper left")
end
ax[2].set_xlabel("θ (degrees)")
ax[1].set_ylim([0,29])
ax[1].set_ylabel(L"\rm  B_c\ (T)")
tight_layout()
subplots_adjust(wspace=0.0)
savefig("test.png",dpi=600,transparent=true)
display(fig)
close(fig)

fig,ax = subplots(figsize=(3,3))
lblstr = ["(-3,-1)","(-2,-2)","(-1,-3)"]
colors = ["tab:blue","tab:red","tab:green","tab:orange","cyan"]
for j in 1:3
    ax.plot(data[:,1],yvals[:,j],"--x",label=lblstr[j],c=colors[j+1])
    ax.legend(loc="upper left")
    ax.set_xlabel("θ (degrees)")
end
ax.set_ylim([0,29])
ax.set_ylabel(L"\rm  B_c\ (T)")
tight_layout()
savefig("test.png",dpi=600,transparent=true)
display(fig)
close(fig)



energies1 = readdlm("120_BM_energies.txt")
energies = readdlm("120_HF_energies.txt")

fig,ax = subplots(1,3,figsize=(5,2.5),sharex=true,sharey=true)
lblstr = ["(-1,-3)","(-2,-2)","(-3,-1)"]
colors = ["tab:blue","tab:red","tab:green","tab:orange","cyan"]
for j in 1:3
    if false # j==3
        ax[j].plot(ϕs[4:end],energies[4:end,4-j]-energies1[4:end,4-j],".-",ms=4,label=lblstr[4-j],c=colors[j+1])
    else
        ax[j].plot(ϕs,energies[:,4-j]-energies1[:,4-j],".-",ms=4,label=lblstr[4-j],c=colors[j+1])
    end
    # ax[j].plot(ϕs,energies1[:,j],"--.",label=lblstr[j],c=colors[j+1])
    ax[j].legend(loc="lower right")
    
end
# ax[1].set_ylim([-0.35,0.06])
# ax[2].set_ylim([-0.3,0.1])
# ax[3].set_ylim([-0.3,0.1])
ax[2].set_xlabel(L"ϕ/ϕ_0")
ax[1].set_ylabel(L"\rm E_{HF}-E_{CHF}\ (meV/u.c.)")
ax[3].set_xlim([0.01,0.53])
ax[3].set_xticks([0.1,0.3,0.5])
tight_layout()
subplots_adjust(wspace=0.0)
savefig("test.png",dpi=600,transparent=false)
display(fig)
close(fig)




energies1 = readdlm("120_BM_energies.txt")
energies = readdlm("120_HF_energies.txt")

fig, ax = subplots(figsize=(3,3))
lblstr = ["(-1,-3)","(-2,-2)","(-3,-1)"]
colors = ["tab:blue","tab:red","tab:green","tab:orange","cyan"]
for j in 1:1
    # ax[j].plot(ϕs,energies[:,j]-energies1[:,j],".-",ms=4,label=lblstr[j],c=colors[j+1])
    ax.plot(ϕs[1:10],energies1[1:10,j],"--.",label=lblstr[j],c=colors[j+1])
    ax.plot(ϕs[1:10],energies[1:10,j],"--.",label=lblstr[j],c="k")
    ax.legend(loc="lower right")
    
end
# ax[1].set_ylim([-0.9,0.08])
ax.set_xlabel(L"ϕ/ϕ_0")
ax.set_ylabel(L"\rm E_{HF}-E_{CHF}\ (meV/u.c.)")
tight_layout()
display(fig)
close(fig)