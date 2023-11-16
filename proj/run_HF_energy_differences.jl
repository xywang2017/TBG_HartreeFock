using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
dir = "/media/xiaoyuw@ad.magnet.fsu.edu/Data/Code/TBG_HartreeFock/"
dir = "/Volumes/Data/Code/TBG_HartreeFock/"
# Info and folder name
# ------------------------------------------------------------------------------ # 
twist_angle = 105
foldername = dir*"zeeman/$(twist_angle)_strain"
# params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)


w0 = "07"
ϕs = sort(unique([p//q for q in 1:12 for p in 1:q]))
ϕs = ϕs[ϕs.<=0.5]
sts = []
for s in -3:3, t in -12:12
    push!(sts,[s,t])
    push!(sts,[-s,-t])
end
sts = unique(sts)

# --------------------------------------- energies --------------------------------------- # 
HF = []
HF_flavor = []
HF_tL2 = []
sts = [[0,-4],[-1,-3],[-2,-2],[-3,-1]]
# sts = [[0,-4]]
for st in sts 
    E_HF = Float64[]
    E_HF_flavor = Float64[]
    E_HF_tL2 = Float64[] 
    s,t = st[1], st[2]
    for ϕ in ϕs 
        p,q = numerator(ϕ), denominator(ϕ)
        νstr = round(Int,1000*(s+t*p/q))
        if abs(s+t*p/q) < 4
            if !(t==-4 )
                metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="bm_cascade_new_init_HF_$(p)_$(q)_nu_$(νstr)")
            else
                metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="init_HF_$(p)_$(q)_nu_$(νstr)")
            end
            if !isempty(metadata)
                push!(E_HF_flavor,load(metadata,"iter_energy")[end])
            else 
                push!(E_HF_flavor,NaN)
            end

            metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="tL_init_HF_$(p)_$(q)_nu_$(νstr)")
            if !isempty(metadata)
                push!(E_HF_tL2,load(metadata,"iter_energy")[end])
            else
                push!(E_HF_tL2,NaN)
            end

            metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="init_HF_$(p)_$(q)_nu_$(νstr)")
            if !isempty(metadata)
                push!(E_HF,load(metadata,"iter_energy")[end])
            else
                push!(E_HF,NaN)
            end
        end
    end
    push!(HF,E_HF)
    push!(HF_flavor,E_HF_flavor)
    push!(HF_tL2,E_HF_tL2)
end

#
ntot = length(HF)
fig, ax = subplots(1,ntot,sharex=true,figsize=(2.8ntot,2.8))

for r in 1:ntot 
    # ax[r].plot(ϕs,HF[r],"ro",label="HF G.S.",ms=3)
    # ax[r].plot(ϕs,HF_flavor[r],"b^",label="CHF",ms=3)
    ax[r].plot(ϕs,HF[r]-HF_flavor[r] ,"-o",ms=5,c="tab:blue")
    # ax[r].plot(ϕs,HF_tL2[r],"g<",label="HF tL2",ms=2)
    if r == 1
        # ax[r].legend(fontsize=8)
    end
    ax[r].set_xlabel(L"ϕ/ϕ_0")
    if r ==1
        ax[r].set_ylabel("HF energy (meV/u.c.)")
        ax[r].set_ylabel(L"\rm E_{G.S.}-E_{no ivc}\ (meV/u.c.)")
    end
    if r ==3
        # ax[r].set_yticks(collect(-22:5:-12)) 
    end
    s, t = sts[r][1], sts[r][2]
    ax[r].set_title("($(s),$(t))")
end

tight_layout()
savefig("test.png",transparent=true,dpi=500)
display(fig)
close(fig)
