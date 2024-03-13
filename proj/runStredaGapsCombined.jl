using PyPlot,JLD2,Printf
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))

dir = "/media/xiaoyuw@ad.magnet.fsu.edu/Data/Code/TBG_HartreeFock/"
dir = "/Volumes/Data/Code/TBG_HartreeFock/"
# dir = "w0_w1_05/"

w0 = "07"
ϕs = sort(unique([p//q for q in 1:12 for p in 1:q]))
ϕs = ϕs[ϕs.<=0.5]
# ϕs = ϕs[ϕs.>0.1]
sts = []
for s in -3:3, t in -12:12
    push!(sts,[s,t])
    push!(sts,[-s,-t])
end
sts = unique(sts)

# -------------------------Streda Line Plot ---------------------------------- # 
fig, ax = subplots(2,3,figsize=(10,6),sharex=true,sharey=true)
twist_angles = [138,132,128,124,120,105]
ax[1,1].set_xlim([-4.1,0.1])
ax[1,1].set_ylim([-0.01,0.59])
for i in eachindex(twist_angles)
    twist_angle = twist_angles[i]
    r, c = (i-1)÷3 + 1, (i-1)%3 + 1
    θ0 = @sprintf "%.2f" twist_angle*0.01
    # ax[r,c].set_title(L"θ=%$(θ0)^\circ")
    ax[r,c].text(-1.15,0.55,L"θ=%$(θ0)^\circ",size=11,c="k")
    # ax[r,c].text(-3.95,0.7,L"θ=%$(θ0)^\circ",size=11,c="k")
    foldername = dir*"zeeman/$(twist_angle)_strain"
    params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
    initParamsWithStrain(params)
    # for lines in -4:0
    #     ax[r,c].axvline(lines,ls=":",c="gray",lw=0.5)
    # end
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
                Δ,_=computegap(metadata;savename="test.png")
                push!(gaps,Δ)
                # if occursin("_tL_",metadata)
                #     println(metadata)
                # end
            end
        end
        ax[r,c].scatter(ns,ones(length(ns))*ϕ,s=gaps.^2/15,c="tab:blue",edgecolor="none")
    end
    if r==1 && c==1 
        ax[r,c].scatter([-3.93],[0.56],s=10^2/15,c="k",edgecolor="none")
        ax[r,c].text(-3.83,0.55,"10 meV",size=11,c="k")
        # ax[r,c].scatter([-1.0],[0.71],s=10^2/10,c="tab:red",edgecolor="none")
        # ax[r,c].text(-0.9,0.7,"10 meV",size=11,c="tab:red")
    end
    ax2 = ax[r,c].twinx()
    mn, mx = ax[r,c].get_ylim()
    ax2.set_ylim(flux_conversion(mn,params), flux_conversion(mx,params))
    ax2.set_yticks(collect(0:5:0.9*flux_conversion(mx,params)))
    if c==3
        ax2.set_ylabel("B (T)",fontsize=13)
    end
end


for c in 1:3 
    ax[2,c].set_xlabel(L"n/n_s",fontsize=13)
end
for r in 1:2
    ax[r,1].set_ylabel(L"ϕ/ϕ_0",fontsize=13)
end
# ax[1,1].set_ylim([0,0.8])
# ax[1,1].set_yticks(collect(0:0.2:0.5))
tight_layout()
savefig(joinpath(fpath,"fig1_strain_color.png"),transparent=true,dpi=600)
display(fig)
close(fig)


# ------------------------------- sublattice polarization analysis versus angle for a given flux --------------------------- 
twist_angles = [105; collect(106:2:138)]
# twist_angles = [105;120;124;128;132;138]
sts = [[0,-4],[-1,-3],[-2,-2],[-3,-1]]
ϕ = 2//5
p,q = numerator(ϕ), denominator(ϕ)
_is_strain = "strain"
# ϕ, s, t =1//8, -1,-3

fig = figure(figsize=(3,2.8))

Pzs_bounds = ComplexF64[]
for i in eachindex(twist_angles)
    twist_angle = twist_angles[i]
    foldername = dir*"zeeman/$(twist_angle)_"*_is_strain
    s, t = 0,-4
    ν = s + t*ϕ
    νstr = round(Int,1000*ν)
    metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="nu_$(νstr)")
    hf = load(metadata,"hf");
    P = reshape(hf.P,8q,8q,:);
    Σz = reshape(hf.Σz0,8q,8q,:);
    nmax = round(Int,(ν+4)/8 *size(P,1)*size(P,3))

    σz_upper = 0.0 
    for ik in 1:size(hf.P,3)
        tmpΣz =  eigvals(Hermitian(Σz[1:(2q),1:(2q),ik]))
        σz_upper += sum(tmpΣz[(q+p+1):(2q)])
    end

    σz_lower = 0.0 
    for ik in 1:size(hf.P,3)
        tmpΣz = diag(Σz[1:(2q),1:(2q),ik])
        σz_lower += sum(tmpΣz[1:(q-p)])
    end
    push!(Pzs_bounds,(σz_upper+1im*σz_lower)/((q-p)*size(P,3)))
end

# plot(twist_angles.*0.01,real(Pzs_bounds),"k:",label="sCI")
# plot(twist_angles.*0.01,imag(Pzs_bounds),"k--",label="HSF")


plot(twist_angles.*0.01,real(Pzs_bounds),"k:")
plot(twist_angles.*0.01,imag(Pzs_bounds),"k--")

for st in sts
    s,t = st[1], st[2]
    Pzs = Float64[]

    for i in eachindex(twist_angles)
        twist_angle = twist_angles[i]
        ν = s + t*ϕ
        νstr = round(Int,1000*ν)
        foldername = dir*"zeeman/$(twist_angle)_"*_is_strain
        # println(s," ",t," ",foldername)
        metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="nu_$(νstr)")
        hf = load(metadata,"hf");
        P = reshape(hf.P,8q,8q,:);
        Σz = reshape(hf.Σz0,8q,8q,:);
        nmax = round(Int,(ν+4)/8 *size(P,1)*size(P,3))

        σzmid = sum([tr(( transpose(view(P,:,:,ik))+0.5I ) * view(Σz,:,:,ik)) /nmax for ik in 1:size(P,3) ] ) 
        push!(Pzs,real(σzmid)) 
    end

    plot(twist_angles.*0.01,Pzs,"-o",ms=2,label="($(s),$(t))")
    # plot(twist_angles.*0.01,Pzs,"-o",ms=2)
end

# ylim([-0.06,0.66])
yticks(collect(0:0.2:0.6))
# legend(loc=[0.55,0.45],fontsize=10)
legend(loc=[0.5,0.4],fontsize=10)
xlabel("θ")
ylabel(L"\rm ⟨σ_zτ_z⟩\ per\ electron")
tight_layout()
savefig("$(_is_strain)_$(p)_$(q).png",dpi=600,transparent=true)
display(fig)
close(fig)

