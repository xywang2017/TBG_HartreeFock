using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
# Hartree Fock related 
params = Params(w1=110,w0=77,vf=2482,dθ=1.05π/180)
initParamsWithStrain(params)

w0 = "07"

ϕs = [1//12;1//10;1//8;1//6;1//5;1//4;2//7;1//3;2//5;3//7;1//2]
sts = unique([[s,t] for s in 0:3 for t in 0:4])
# ------------------------------------------ # 
for st in sts 
    s,t = st[1], st[2]
    # s,t = 0,0
    metadatas = String[]
    for ϕ in ϕs 
        p,q = numerator(ϕ), denominator(ϕ)
        νstr = round(Int,1000*(s+t*p/q))
        if s+t*p/q < 4
            metadata = joinpath(fpath,"feldman/B/data_w07/_$(p)_$(q)/1_random_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
            E = load(metadata,"iter_energy")[end]
            for flag in ["flavor","random","chern"], seed in 1:6 
                metadata0 = joinpath(fpath,"feldman/B/data_w$(w0)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
                if isfile(metadata0)
                    E0 = load(metadata0,"iter_energy")[end]
                    if E0<=E 
                        E, metadata = E0, metadata0 
                    end
                end
            end
            println("Iter err: ",load(metadata,"iter_err")[end])
            # println(metadata)
            push!(metadatas,metadata)
        end
    end
    Δs= plot_spectra_collective(metadatas;savename="spectrum_s$(s)_t$(t).png",titlestr="(s,t)=($(s),$(t))");
end


## -------------------------------- plot all gaps ------------------------------ #
ϕs = [1//12;1//10;1//8;1//6;1//5;1//4;2//7;1//3;2//5;3//7;1//2]
# ϕs = [1//8;1//4;2//7;3//7;1//2]
cs = ["r";"g";"b";"c";"m";"darkviolet";"tab:blue";
        "magenta";"peru";"tab:purple";"tab:olive";"deepskyblue";"seagreen";"gray"]
fig = figure(figsize=(6,4))
for lines in [0.0,1.0,2.0,3.0,4.0]
    axvline(lines,ls=":",c="gray")
end
for ϕ in ϕs 
    p,q = numerator(ϕ), denominator(ϕ)
    gaps = Float64[]
    fillings = Float64[]
    for st in sts 
        s,t = st[1], st[2]
        νstr = round(Int,1000*(s+t*p/q))
        if s+t*p/q < 4
            push!(fillings,s+t*p/q)
            metadata = joinpath(fpath,"feldman/B/data_w07/_$(p)_$(q)/1_random_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
            E = load(metadata,"iter_energy")[end]
            for flag in ["flavor","random","chern"], seed in 1:6 
                metadata0 = joinpath(fpath,"feldman/B/data_w07/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
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
    idx_sort = sortperm(fillings)
    fillings, gaps = fillings[idx_sort],gaps[idx_sort]
    idx = unique(z -> fillings[z], 1:length(fillings))
    fillings,gaps = fillings[idx], gaps[idx]
    scatter(fillings,ones(length(fillings))*ϕ,s=gaps.^2/5,c="gray")
end
xlim([-0.3,4.3])
ylim([0.0,0.55])
xlabel("ν")
ylabel(L"ϕ/ϕ_0")
tight_layout()
savefig("streda_line.pdf",transparent=true)
display(fig)
close(fig)