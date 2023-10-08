using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
dir = "/media/xiaoyuw@ad.magnet.fsu.edu/Data/Code/TBG_HartreeFock/"
dir = "/Volumes/Data/Code/TBG_HartreeFock/"
# Info and folder name
# ------------------------------------------------------------------------------ # 
twist_angle = 120
foldername = dir*"zeeman/$(twist_angle)_nostrain"
params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)


w0 = "07"

ϕs = sort(unique([p//q for q in 1:12 for p in 1:q]))
ϕs = ϕs[ϕs.<=0.5]
ϕs = ϕs[ϕs .>0.1]
sts = [[0,-4],[-1,-3],[-2,-2],[-3,-1]]
# -------------------------Streda Line Plot ---------------------------------- # 
cs = ["r";"g";"b";"c";"m";"darkviolet";"tab:blue";
        "magenta";"peru";"tab:purple";"tab:olive";"deepskyblue";"seagreen";"gray"]
fig = figure(figsize=(3,2.2))
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
    colors = []
    for ν in fillings 
        νstr = round(Int,1000*ν)
        metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="nu_$(νstr)")
        purple = 0.5*[1,0,1]
        red = [1,0,0] * 0.5
        green = [0,1,0] *0.8
        if !isempty(metadata)
            push!(ns,ν)
            Δ,Pz=computegap(metadata;savename="test.png")
            push!(gaps,Δ)
            # if (ϕ<1//3) && abs(ν-(-3-ϕ))<1e-2
            #     push!(colors,green)
            # elseif (ϕ<3//10) && abs(ν-(-2-2ϕ))<1e-2
            #     push!(colors,green)
            # elseif (ϕ<2//5) && abs(ν-(-1-3ϕ))<1e-2
            #     push!(colors,green)
            # else
            #     push!(colors,red)
            # end
            if (ϕ<1//6) && abs(ν-(-3-ϕ))<1e-2
                push!(colors,purple)
            elseif (ϕ<1//10) && abs(ν-(-2-2ϕ))<1e-2
                push!(colors,purple)
            elseif (ϕ<1//10) && abs(ν-(-1-3ϕ))<1e-2
                push!(colors,purple)
            else
                push!(colors,red)
            end
        end
    end
    scatter(ns,ones(length(ns))*ϕ,s=gaps.^2/10,c=colors,edgecolor="none")
end
xlim([-4.3,0.3])
ylim([0.0,0.55])
xlabel(L"n/n_s")
ylabel(L"ϕ/ϕ_0")
tight_layout()
savefig(joinpath(fpath,"streda_line_$(twist_angle).png"),transparent=true,dpi=600)
display(fig)
close(fig)
