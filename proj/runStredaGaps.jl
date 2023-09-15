using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
# Info and folder name
# ------------------------------------------------------------------------------ # 
twist_angle = 128
foldername = "zeeman/$(twist_angle)_strain"
params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)


w0 = "07"

# ϕs = [1//8;1//6;1//5;1//4;2//7;1//3;2//5;3//7;1//2]
ϕs = sort(unique([p//q for q in 1:12 for p in 1:q]))
ϕs = ϕs[ϕs.<=0.5]
ϕs = ϕs[2:end]
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
    fillings = Float64[]
    for st in sts 
        s,t = st[1], st[2]
        νstr = round(Int,1000*(s+t*p/q))
        if abs(s+t*p/q) < 4 && (s+t*p/q) <=0
            metadata = joinpath(fpath,"$(foldername)/_$(p)_$(q)/1_random_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
            if !isfile(metadata)
                metadata = joinpath(fpath,"$(foldername)/_$(p)_$(q)/1_flavor_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
            end
            if isfile(metadata)
                push!(fillings,s+t*p/q)
                E = load(metadata,"iter_energy")[end]
                for flag in ["flavor","random","chern","bm","strong","bm_cascade"], seed in 1:10 
                    metadata0 = joinpath(fpath,"$(foldername)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
                    if isfile(metadata0)
                        E0 = load(metadata0,"iter_energy")[end]
                        if E0<=E 
                            E, metadata = E0, metadata0 
                        end
                    end
                end 
                Δ=computegap(metadata;savename="test.pdf")
                push!(gaps,Δ)
            end
        end
    end
    fillings = round.(fillings,digits=8)
    idx_sort = sortperm(fillings)
    fillings, gaps = fillings[idx_sort],gaps[idx_sort]
    idx = unique(z -> fillings[z], 1:length(fillings))
    fillings,gaps = fillings[idx], gaps[idx]
    
    scatter(fillings,ones(length(fillings))*ϕ,s=gaps.^2/10,c="k",edgecolor="none")
    # for t in 1:4 
    #     s = -4 
    #     νtest = round(s + (p/q)*t,digits=8)
    #     idx = fillings .==νtest
    #     scatter(fillings[idx],ones(length(fillings[idx]))*ϕ,s=gaps[idx].^2/5,c="b",edgecolor="none")
    # end
end
xlim([-4.3,0.3])
ylim([0.0,0.55])
xlabel(L"n/n_s")
ylabel(L"ϕ/ϕ_0")
tight_layout()
savefig(joinpath(fpath,"$(foldername)/streda_line.png"),transparent=false,dpi=500)
display(fig)
close(fig)


# -----------------------------Hofstadter spectrum plot ---------------------------- # 
Δss = []
sts = [[-3,-1]]
for st in sts 
    s,t = st[1], st[2]
    metadatas = String[]
    for ϕ in ϕs 
        p,q = numerator(ϕ), denominator(ϕ)
        νstr = round(Int,1000*(s+t*p/q))
        if abs(s+t*p/q) < 4
            metadata = joinpath(fpath,"$(foldername)/_$(p)_$(q)/1_random_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
            if !isfile(metadata)
                metadata = joinpath(fpath,"$(foldername)/_$(p)_$(q)/1_flavor_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
            end
            if isfile(metadata)
                E = load(metadata,"iter_energy")[end]
                for flag in ["flavor","random","chern","bm","strong","bm_cascade"], seed in 1:10
                    metadata0 = joinpath(fpath,"$(foldername)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
                    if isfile(metadata0)
                        E0 = load(metadata0,"iter_energy")[end]
                        if E0<=E 
                            E, metadata = E0, metadata0 
                        end
                    end
                end
                if load(metadata,"iter_err")[end] > 1e-6
                    println("s= ",s," t=",t," p=",p," q=",q," Iter err: ",load(metadata,"iter_err")[end])
                end
                # println(metadata)
                push!(metadatas,metadata)
            end
        end
    end
    Δs= plot_spectra_collective(metadatas;savename="spectrum_s$(s)_t$(t).png",titlestr="(s,t)=($(s),$(t))");
    push!(Δss,Δs)
end

fig = figure(figsize=(3,2.5))
sts = [[0,-4],[-1,-3],[-2,-2],[-3,-1]]
for i in eachindex(Δss)
    s, t = sts[i][1], sts[i][2]
    plot(ϕs,Δss[i],"o-",ms=4,markeredgecolor="none",label="($(s),$(t))")
end
legend(fontsize=8)
xlabel(L"ϕ/ϕ_0")
ylabel("Δ (meV)")
# axvline(3/8)
# axvline(2/7)
# axvline(3/7)
# axvline(4/11)
xticks([0.1,0.2,0.3,0.4,0.5])
xlim([0,0.55])
# ylim([0,1.2*maximum(Δs)])
tight_layout()
savefig("tmp_128.png",transparent=true,dpi=500)
display(fig)
close(fig)

# ---------------- non interacting hofstadter spectrum weighted with a given Streda line density matrix
# plot spectrum 
function plot_LL_spectrum(s::Int,t::Int,params::Params)
    foldername0 = "NonInt/$(twist_angle)_strain"
    fig,ax = subplots(4,1,sharex=true,sharey=true,figsize=(2.5,5))
    strs = ["K","Kprime"]
    for ϕ in ϕs
        p, q = numerator(ϕ), denominator(ϕ)
        νstr = round(Int,1000*(s+t*p/q))
        metadata = joinpath(fpath,"$(foldername)/_$(p)_$(q)/1_random_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
        if !isfile(metadata)
            metadata = joinpath(fpath,"$(foldername)/_$(p)_$(q)/1_flavor_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
        end
        if isfile(metadata)
            E = load(metadata,"iter_energy")[end]
            for flag in ["flavor","random","chern","bm","strong","bm_cascade"], seed in 1:10
                metadata0 = joinpath(fpath,"$(foldername)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
                if isfile(metadata0)
                    E0 = load(metadata0,"iter_energy")[end]
                    if E0<=E 
                        E, metadata = E0, metadata0 
                    end
                end
            end
        end
        hf = load(metadata,"hf");
        P = reshape(hf.P,2hf.q,hf.nη,hf.ns,2hf.q,hf.nη,hf.ns,hf.q,hf.nq,hf.nq)
        weights = zeros(Float64,2hf.q,hf.nq,hf.nq)
        titlestr = ["K↑" "K↓";"K'↑" "K'↓"]
        for iη in 1:2, is in 1:2
            r = (is-1)*2+iη
            str = strs[iη]
            fname0 = joinpath(fpath,"$(foldername0)/_$(p)_$(q)_$(str)_metadata.jld2")
            energies = load(fname0,"E");
            for i2 in 1:hf.nq, i1 in 1:hf.nq, iq in 1:hf.q
                weights[:,i1,i2] = real(diag(P[:,iη,is,:,iη,is,iq,i1,i2])) .+0.5
            end
            energies = reshape(energies,2q,:)
            ax[r].scatter(ones(length(energies[1:(q-p),:][:]))*ϕ,energies[1:(q-p),:][:].+(3-2is)*p/q*ZeemanUnit(params),marker="o",s=3,edgecolor="none",c="b",vmin=0,vmax=1,cmap="coolwarm")
            ax[r].scatter(ones(length(energies[(q-p+1):(2q),:][:]))*ϕ,energies[(q-p+1):(2q),:][:].+(3-2is)*p/q*ZeemanUnit(params),marker="o",s=3,edgecolor="none",c="gray",vmin=0,vmax=1,cmap="coolwarm")
            # ax[iη,is].set_title(titlestr[iη,is])
            ax[r].text(0.02,0.0,titlestr[iη,is],rotation="vertical")
        end
    end
    for j in 1:4
        ax[j].set_ylabel("E (meV)")
        if j==4
            ax[j].set_xlabel(L"ϕ/ϕ_0")
        end
    end
    ax[1,1].set_xlim([0,0.55])
    tight_layout() 
    subplots_adjust(hspace=0,wspace=0)
    savefig("tmp_$(s)_$(t).png",dpi=500)
    display(fig)
    close(fig)
    return nothing
end

plot_LL_spectrum(-1,-3,params)

