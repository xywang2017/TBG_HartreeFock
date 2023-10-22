using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
# Info and folder name
# ------------------------------------------------------------------------------ # 
dir = "/media/xiaoyuw@ad.magnet.fsu.edu/Data/Code/TBG_HartreeFock/"
# dir = "/Volumes/Xiaoyu/Code/TBG_HartreeFock/"
# dir = ""
twist_angle = 138
foldername = dir*"zeeman/$(twist_angle)_strain"
params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)

w0 = "07"

ϕ = 1//8
sts = []
for s in -3:3, t in -12:12
    push!(sts,[s,t])
    push!(sts,[-s,-t])
end
sts = unique(sts)
νs = sort(unique([sts[i][1]+sts[i][2]*ϕ*1.0 for i in eachindex(sts)]))
νs = νs[abs.(νs) .<=4]
# νs = νs[νs .<=1e-3]
#
fig = figure(figsize=(3,3))
ϵs = collect(range(-95,95,400))
γ = 1
dos = zeros(Float64,length(ϵs),length(νs))
p,q = numerator(ϕ), denominator(ϕ)
for j in eachindex(νs)
    ν = νs[j]
    νstr = round(Int,1000*ν)
    metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="_init_HF_$(p)_$(q)_nu_$(νstr)",_printinfo=false)
    if isfile(metadata)
        energies = load(metadata,"hf").ϵk[:]
        μ =load(metadata,"hf").μ
        # plot(energies,ones(length(energies))*ν,"ko",ms=2,markeredgecolor="none")
        for i in eachindex(ϵs)
            dos[i,j] = 1/π * sum(γ./(γ^2 .+ (ϵs[i].-(energies.-μ)).^2)) / length(energies)
        end
    end
end
pl = pcolor(ϵs,νs,dos' ./ maximum(dos),cmap="bwr",vmin=0,vmax=1)
# for sts in [[-1,-3],[-2,-2],[-3,-1],[1,3],[2,2],[3,1],[0,4],[0,-4]]
#     s,t = sts[1],sts[2]
#     axhline(s+t*ϕ,ls=":",c="gray")
# end

# axhline(-2.875,ls=":",c="k")
colorbar(pl,shrink=0.8)
ylim([-4,4])
yticks(collect(-3:3))
ylabel(L"n/n_s")
xlabel("E (meV)")
tight_layout()
display(fig)
savefig("$(twist_angle)_strain_spectra.png",dpi=500,transparent=true)
# savefig("main_lines.png",dpi=500,transparent=true)
close(fig)