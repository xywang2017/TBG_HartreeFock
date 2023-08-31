using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
# Info and folder name
# ------------------------------------------------------------------------------ # 
twist_angle = 105
foldername = "zeeman/$(twist_angle)_strain"
params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)


w0 = "07"

ϕs = [1//8;1//6;1//5;1//4;2//7;1//3;2//5;3//7;1//2]
# ϕs = [1//12;1//10;1//8;1//6]
sts = []
for s in -3:3, t in -12:12
    push!(sts,[s,t])
    push!(sts,[-s,-t])
end
sts = unique(sts)

# -------------------------Streda Line Plot ---------------------------------- # 
cs = ["r";"g";"b";"c";"m";"darkviolet";"tab:blue";
        "magenta";"peru";"tab:purple";"tab:olive";"deepskyblue";"seagreen";"gray"]
fig = figure(figsize=(5,4))
for lines in -4:4
    axvline(lines,ls=":",c="gray")
end
for ϕ in ϕs 
    p,q = numerator(ϕ), denominator(ϕ)
    gaps = Float64[]
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
            end
        end
    end
    # fillings = round.(fillings,digits=8)
    idx_sort = sortperm(fillings)
    fillings, gaps = fillings[idx_sort],gaps[idx_sort]
    idx = unique(z -> fillings[z], 1:length(fillings))
    fillings,gaps = fillings[idx], gaps[idx]
    
    scatter(fillings,ones(length(fillings))*ϕ,s=gaps.^2/5,c="k",edgecolor="none")
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
sts = [[-1.5,-2]]
for st in sts 
    s,t = st[1], st[2]
    metadatas = String[]
    for ϕ in ϕs 
        p,q = numerator(ϕ), denominator(ϕ)
        νstr = round(Int,1000*(s+t*p/q))
        if abs(s+t*p/q) < 4
            metadata = joinpath(fpath,"$(foldername)/B/data_w07/_$(p)_$(q)/1_random_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
            if !isfile(metadata)
                metadata = joinpath(fpath,"$(foldername)/B/data_w07/_$(p)_$(q)/1_flavor_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
            end
            if isfile(metadata)
                E = load(metadata,"iter_energy")[end]
                for flag in ["flavor","random","chern","bm","strong","bm_cascade"], seed in 1:10
                    metadata0 = joinpath(fpath,"$(foldername)/B/data_w$(w0)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
                    if isfile(metadata0)
                        E0 = load(metadata0,"iter_energy")[end]
                        if E0<=E 
                            E, metadata = E0, metadata0 
                        end
                    end
                end
                if load(metadata,"iter_err")[end] > 1e-6
                    println("s= ",s," t=",t," p=",p," q=",q," Iter err: ",load(metadata,"iter_energy")[end])
                end
                println(metadata)
                push!(metadatas,metadata)
            end
        end
    end
    Δs= plot_spectra_collective(metadatas;savename="spectrum_s$(s)_t$(t).png",titlestr="(s,t)=($(s),$(t))");

    fig = figure(figsize=(3,2))
    plot(ϕs,Δs,"b^-",ms=5,markeredgecolor="none",label="($(s),$(t))")
    legend()
    xlabel(L"ϕ/ϕ_0")
    ylabel("Δ (meV)")
    xlim([0,0.53])
    ylim([0,1.2*maximum(Δs)])
    tight_layout()
    savefig("tmp.png",transparent=true,dpi=500)
    display(fig)
    close(fig)
end
