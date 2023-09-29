using PyPlot,JLD2,Printf
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))


w0 = "07"
ϕs = sort(unique([p//q for q in 1:12 for p in 1:q]))
ϕs = ϕs[ϕs.<=0.5]
ϕs = ϕs[ϕs.>0.1]
sts = []
for s in -3:3, t in -12:12
    push!(sts,[s,t])
    push!(sts,[-s,-t])
end
sts = unique(sts)

# -------------------------Streda Line Plot ---------------------------------- # 
fig, ax = subplots(2,3,figsize=(10,6),sharex=true,sharey=true)
twist_angles = [138,132,128,124,120,105]
for i in eachindex(twist_angles)
    twist_angle = twist_angles[i]
    r, c = (i-1)÷3 + 1, (i-1)%3 + 1
    θ0 = @sprintf "%.2f" twist_angle*0.01
    # ax[r,c].set_title(L"θ=%$(θ0)^\circ")
    ax[r,c].text(-3.95,0.55,L"θ=%$(θ0)^\circ",size=11,c="k")
    foldername = "zeeman/$(twist_angle)_strain"
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
            end
        end
        ax[r,c].scatter(ns,ones(length(ns))*ϕ,s=gaps.^2/10,c="k",edgecolor="none")
    end
    if r==1 && c==1 
        ax[r,c].scatter([-1.0],[0.56],s=10^2/10,c="tab:red",edgecolor="none")
        ax[r,c].text(-0.9,0.55,"10 meV",size=11,c="tab:red")
    end
end

ax[1,1].set_xlim([-4.1,0.1])
ax[1,1].set_ylim([0.05,0.59])

for c in 1:3 
    ax[2,c].set_xlabel(L"n/n_s")
end
for r in 1:2
    ax[r,1].set_ylabel(L"ϕ/ϕ_0")
end
tight_layout()
savefig(joinpath(fpath,"fig1_strain.png"),transparent=false,dpi=600)
display(fig)
close(fig)



# ------------------------------- sublattice polarization analysis versus angle for a given flux --------------------------- 
twist_angles = [138,132,128,124,120,105]
ϕ, s, t = 2//5, -3, -1
Pzs_strain = Float64[]
Pzs_strain_bounds = ComplexF64[]
Pzs_nostrain = Float64[]
Pzs_nostrain_bounds = ComplexF64[]

for i in eachindex(twist_angles)
    twist_angle = twist_angles[i]
    p,q = numerator(ϕ), denominator(ϕ)
    ν = s + t*ϕ
    νstr = round(Int,1000*ν)

    foldername0 = "NonInt/$(twist_angle)_nostrain/_$(p)_$(q)_K_metadata.jld2"
    PΣz = reshape(load(foldername0,"PΣz"),2q,2q,:);
    σzlower = 0.0 
    for i in 1:(q-p)
        σzlower += sum(PΣz[i,i,:])
    end
    σzlower /= ( size(PΣz,3)*(q-p) )
    σzupper = 0.0
    for i in 1:size(PΣz,3)
        vals = eigvals(Hermitian(PΣz[:,:,i]))
        σzupper += sum(vals[(q+p+1):(2q)])
    end
    σzupper /=( size(PΣz,3)*(q-p) )
    push!(Pzs_nostrain_bounds,σzlower+1im*σzupper)

    foldername = "zeeman/$(twist_angle)_nostrain"
    metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="nu_$(νstr)")

    if !isempty(metadata)
        Δ,Pz=computegap(metadata;savename="test.png")
        push!(Pzs_nostrain,Pz)
    end


    foldername0 = "NonInt/$(twist_angle)_strain/_$(p)_$(q)_K_metadata.jld2"
    PΣz = reshape(load(foldername0,"PΣz"),2q,2q,:);
    σzlower = 0.0 
    for i in 1:(q-p)
        σzlower += sum(PΣz[i,i,:])
    end
    σzlower /= ( size(PΣz,3)*(q-p) )
    σzupper = 0.0
    for i in 1:size(PΣz,3)
        vals = eigvals(Hermitian(PΣz[:,:,i]))
        σzupper += sum(vals[(q+p+1):(2q)])
    end
    σzupper /=( size(PΣz,3)*(q-p) )
    push!(Pzs_strain_bounds,σzlower+1im*σzupper)
    

    foldername = "zeeman/$(twist_angle)_strain"
    metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="nu_$(νstr)")
    if !isempty(metadata)
        Δ,Pz=computegap(metadata;savename="test.png")
        push!(Pzs_strain,Pz)
    end
end

fig = figure(figsize=(4,3))
plot(twist_angles.*0.01,real(Pzs_nostrain_bounds),":",c="r")
plot(twist_angles.*0.01,imag(Pzs_nostrain_bounds),":",c="r")
plot(twist_angles.*0.01,Pzs_nostrain,"r-v",label="($s,$t) no strain")

plot(twist_angles.*0.01,real(Pzs_strain_bounds),":",c="b")
plot(twist_angles.*0.01,imag(Pzs_strain_bounds),":",c="b")
plot(twist_angles.*0.01,Pzs_strain,"b-^",label="($s,$t) strain")
legend(loc="upper right")
xlabel("θ")
ylabel(L"\rm ⟨σ_zτ_z⟩")
tight_layout()
savefig("Pz_s_$(s)_t_$(t).png",dpi=600,transparent=true)
display(fig)
close(fig)



foldername0 = "NonInt/105_strain/_2_5_K_metadata.jld2"
PΣz = reshape(load(foldername0,"PΣz"),2q,2q,:);
σzlower = 0.0 
for i in 1:(q-p)
    σzlower += sum(PΣz[i,i,:])
end
σzlower /= ( size(PΣz,3)*(q-p) )
σzupper = 0.0
for i in 1:size(PΣz,3)
    vals = eigvals(Hermitian(PΣz[:,:,i]))
    σzupper += sum(vals[(q+p+1):(2q)])
    if i==4
        fig = figure()
        plot(vals,"bo")
        display(fig)
        close(fig) 
    end
end