using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
# Info and folder name
# ------------------------------------------------------------------------------ # 
twist_angle = 120
foldername = "zeeman/$(twist_angle)_nostrain"
params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)


w0 = "07"

# ϕs = [1//8;1//6;1//5;1//4;2//7;1//3;2//5;3//7;1//2]
ϕs = sort(unique([p//q for q in 1:12 for p in 1:q]))
ϕs = ϕs[ϕs.<=0.5]
ϕs = ϕs[4:end]
sts = []
for s in -3:3, t in -12:12
    push!(sts,[s,t])
    push!(sts,[-s,-t])
end
sts = unique(sts)
# -------------------------Streda Line Plot ---------------------------------- # 
cs = ["r";"g";"b";"c";"m";"darkviolet";"tab:blue";
        "magenta";"peru";"tab:purple";"tab:olive";"deepskyblue";"seagreen";"gray"]
fig = figure(figsize=(3,3))
for lines in -4:4
    axvline(lines,ls=":",c="gray",lw=0.5)
end
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
            Δ=computegap(metadata;savename="test.png")
            push!(gaps,Δ)
        end
    end
    scatter(ns,ones(length(ns))*ϕ,s=gaps.^2/10,c="k",edgecolor="none")
end
xlim([-4.3,0.3])
ylim([0.0,0.55])
xlabel(L"n/n_s")
ylabel(L"ϕ/ϕ_0")
tight_layout()
savefig(joinpath(fpath,"$(foldername)/streda_line.png"),transparent=false,dpi=600)
display(fig)
close(fig)


# -----------------------------Hofstadter spectrum plot ---------------------------- # 
Δss = []
# sts = [[0,0],[0,-3],[0,-2],[0,-1]]
sts = [[-3,-1]]
for st in sts 
    s,t = st[1], st[2]
    metadatas = String[]
    for ϕ in ϕs 
        p,q = numerator(ϕ), denominator(ϕ)
        νstr = round(Int,1000*(s+t*p/q))
        if abs(s+t*p/q) < 4
            metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="nu_$(νstr)")
            if !isempty(metadata)
                push!(metadatas,metadata)
            end
        end
    end
    Δs= plot_spectra_collective(metadatas;savename="spectrum_s$(s)_t$(t).png",titlestr="(s,t)=($(s),$(t))");
    push!(Δss,Δs)
end

fig = figure(figsize=(6,4))
# sts = [[0,0],[0,-3],[0,-2],[0,-1]]
# ϕc = [-1,2//5,1//3,1//3] # 132 degrees critical field
# ϕc = [-1,1//6,1//11,1//11] # 132 degrees critical field
colors = ["tab:blue","tab:red","tab:green","tab:orange","cyan"]
for i in eachindex(Δss)
    s, t = sts[i][1], sts[i][2]
    plot(ϕs,Δss[i],"-o",c=colors[i],ms=4,markeredgecolor="none",label="($(s),$(t))")
end
legend(fontsize=6,loc="upper left")
xlabel(L"ϕ/ϕ_0")
ylabel("Δ (meV)")
# axvline(3/8)
# axvline(2/7)
# axvline(3/7)
# axvline(4/11)
# xticks([0.1,0.2,0.3,0.4,0.5])
xlim([-0.1,0.55])
ylim([0,24])
tight_layout()
savefig("tmp_128.png",transparent=false,dpi=600)
display(fig)
close(fig)

# ---------------- non interacting hofstadter spectrum weighted with a given Streda line density matrix
# plot spectrum 
function plot_LL_spectrum(params::Params;angle::Int=120)
    foldername0 = "NonInt/$(angle)_strain"
    fig,ax = subplots(figsize=(2.8,2.5))
    # fig,ax = subplots(figsize=(6,4))
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
            push!(tmp,energies[:,1])
            # pl = ax.scatter(ones(length(energies[:]))*ϕ,energies[:],marker="o",s=3,edgecolor="none",c=colors[iη])
            pl = ax.scatter(ones(length(energies[:]))*ϕ,energies[:],marker="o",s=2,edgecolor="none",c=weights[:]*(3-2iη),vmin=-1,vmax=1,cmap="coolwarm")
        end
        
    end
    colorbar(pl,fraction=0.06, pad=0.05)
    ax.set_ylabel("E (meV)")
    ax.set_xlabel(L"ϕ/ϕ_0")
    ax.set_xlim([0,0.55])
    # ax.set_ylim([-39,39])
    tight_layout()
    # savefig("tmp_$(angle).png",dpi=600,transparent=true)
    display(fig)
    close(fig)
    return nothing
end

plot_LL_spectrum(params,angle=120)

