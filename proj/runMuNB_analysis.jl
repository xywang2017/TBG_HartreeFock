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

ϕs = [1//8;1//6;1//5;1//4;2//7;1//3;2//5;3//7;1//2]
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
fig,ax = subplots(figsize=(5,4))
# for lines in -4:4
#     axvline(lines,ls=":",c="gray")
# end
for ϕ in ϕs 
    p,q = numerator(ϕ), denominator(ϕ)
    gaps = Float64[]
    μs = Float64[] 
    fillings = Float64[]
    for st in sts 
        s,t = st[1], st[2]
        νstr = round(Int,1000*(s+t*p/q))
        if abs(s+t*p/q) < 4 && (s+t*p/q) <=0
            metadata = joinpath(fpath,"$(foldername)/B/data_w07/_$(p)_$(q)/1_random_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
            if !isfile(metadata)
                metadata = joinpath(fpath,"$(foldername)/B/data_w07/_$(p)_$(q)/1_flavor_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
            end
            if isfile(metadata)
                push!(fillings,s+t*p/q)
                E = load(metadata,"iter_energy")[end]
                for flag in ["flavor","random","chern","bm","strong","bm_cascade"], seed in 1:10 
                    metadata0 = joinpath(fpath,"$(foldername)/B/data_w07/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
                    if isfile(metadata0)
                        E0 = load(metadata0,"iter_energy")[end]
                        if E0<=E 
                            E, metadata = E0, metadata0 
                        end
                    end
                end 
                Δ=computegap(metadata;savename="test.pdf")
                push!(gaps,Δ)
                push!(μs,load(metadata,"hf").μ)
            end
        end
    end
    fillings = round.(fillings,digits=8)
    idx_sort = sortperm(fillings)
    fillings, gaps, μs = fillings[idx_sort],gaps[idx_sort], μs[idx_sort]
    idx = unique(z -> fillings[z], 1:length(fillings))
    fillings,gaps,μs = fillings[idx], gaps[idx],μs[idx]
    ns = (fillings .+4)./4
    ax.scatter(ones(length(μs))*ϕ,μs,c=ns,vmin=0,vmax=1,cmap="coolwarm",s=3)
    # for i in eachindex(μs)
    #     c_ns = 1 - ns[i]
    #     # ax.plot([ϕ,ϕ],[μs[i]-gaps[i]/2,μs[i]+gaps[i]/2],"o",c=c_ns,ms=4,markeredgecolor="none")
    #     if gaps[i] < 100.0  # only draw gaples states
    #         ax.plot([ϕ],[μs[i]],"o",c=[c_ns,0.2,0.2],ms=4,markeredgecolor="none")
    #     else # draw states above and below the gap
    #         ax.plot([ϕ,ϕ],[μs[i]-gaps[i]/2,μs[i]+gaps[i]/2],"o",c=c_ns,ms=4,markeredgecolor="none")
    #     end
    # end
    # ax.scatter(ones(length(fillings))*ϕ,μs,s=gaps.^2/4,c="k",edgecolor="none")
end
ax.set_xlim([0.0,0.55])
ax.set_ylabel("E (meV)")
ax.set_xlabel(L"ϕ/ϕ_0")
tight_layout()
# savefig(joinpath(fpath,"$(foldername)/streda_line.png"),transparent=false,dpi=500)
display(fig)
close(fig)
