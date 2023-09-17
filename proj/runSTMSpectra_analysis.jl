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

ϕ = 5//12
sts = []
for s in -3:3, t in -12:12
    push!(sts,[s,t])
    push!(sts,[-s,-t])
end
sts = unique(sts)
νs = sort(unique([sts[i][1]+sts[i][2]*ϕ*1.0 for i in eachindex(sts)]))
νs = νs[abs.(νs) .<=4]
νs = νs[νs .<=1e-3]
#
fig = figure(figsize=(3,2))
ϵs = collect(range(-90,90,400))
γ = 2
dos = zeros(Float64,length(ϵs),length(νs))
p,q = numerator(ϕ), denominator(ϕ)
for j in eachindex(νs)
    ν = νs[j]
    νstr = round(Int,1000*ν)
    if ν < 5
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
            energies = load(metadata,"hf").ϵk[:]
            μ =load(metadata,"hf").μ
            # plot(energies,ones(length(energies))*ν,"ko",ms=2,markeredgecolor="none")
            for i in eachindex(ϵs)
                dos[i,j] = 1/π * sum(γ./(γ^2 .+ (ϵs[i].-(energies.-μ)).^2))
            end
        end
    end
end
pl = pcolor(ϵs,νs,dos',cmap="bwr")
# for sts in [[-1,-3],[-2,-2],[-3,-1],]
#     s,t = sts[1],sts[2]
#     axhline(s+t*ϕ,ls=":",c="k")
# end

# axvline(0,ls=":",c="k")
# colorbar(pl)
ylabel(L"n/n_s")
xlabel("E (meV)")
tight_layout()
display(fig)
savefig("$(twist_angle)_strain_spectra.png",dpi=500)
close(fig)