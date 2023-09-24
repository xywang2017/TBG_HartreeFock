using PyPlot 

## plot Hartree Fock spectra
function plot_spectra(metadata::String;savename::String="tmp.pdf")
    hf = load(metadata,"hf");
    νF = 8*round(Int,(hf.ν+4)/8*length(hf.ϵk)) / length(hf.ϵk)-4
    ϵk = hf.ϵk
    σzτz = hf.σzτz
    params = hf.params 

    ee = 1.6e-19
    ϵϵ = 8.8541878128e−12	
    aa = 2.46e-10
    ϵr = 15.0
    Vcoulomb = ee/(4π*ϵϵ*ϵr* abs(params.a1)*aa) * 1e3
    println("Coulomb scale: ",Vcoulomb)
    fig = figure(figsize=(3,3))
    idx = sortperm(ϵk[:])
    ϵsorted = ϵk[idx] #./Vcoulomb
    chern = σzτz[idx]

    # ϵsorted = ϵsorted[chern .>0] #./Vcoulomb
    # chern = chern[chern .>0]

    pl=scatter(ones(length(ϵsorted))*hf.p/hf.q,ϵsorted,c=chern,cmap="coolwarm",s=6,vmin=-1,vmax=1)
    colorbar(pl)
    ν = eachindex(ϵsorted) ./ length(ϵsorted)
    i = 1
    while (νF+4)/8 > ν[i]
        i += 1
    end
    if i<length(chern)
        # i = i-1
        ϵF = (ϵsorted[i+1] + ϵsorted[i])/2 
        Δ = (ϵsorted[i+1] - ϵsorted[i]) 
    else
        ϵF = ϵsorted[end]
        Δ = 0 
    end
    # Δ = (ϵsorted[i] - ϵsorted[i-1]) 
    println("Gap size: ", Δ)
    axhline(ϵF,ls=":",c="gray")
    ylabel("E (meV)")
    xlabel(L"ϕ/ϕ_0")
    # legend()
    # ylim([-0.4,0.8])
    # ylim([-25,35])
    tight_layout()
    # savefig(savename,transparent=true)
    display(fig)
    close(fig)
    # println("Sublattice polarization operator is: ",sum(chern[ϵsorted.<ϵF])/(size(ϵk,2)*size(ϵk,1))*8)
    # println("Total energy: ",(0.25*sum(ϵsorted[ϵsorted.<ϵF])-0.25*sum(ϵsorted[ϵsorted.>=ϵF]))/size(hf.ϵk,2)/size(hf.ϵk,1)*8)
    # return sum(chern[ϵsorted.<ϵF])/(size(ϵk,2)*size(ϵk,1))*8
    return Δ
end


## Compute spectral gap from HF 
function computegap(metadata::String;savename::String="tmp.pdf")
    hf = load(metadata,"hf");
    νF = 8*round(Int,(hf.ν+4)/8*length(hf.ϵk)) / length(hf.ϵk)-4
    ϵk = hf.ϵk
    σzτz = hf.σzτz
    params = hf.params 
    idx = sortperm(ϵk[:])
    ϵsorted = ϵk[idx] #./Vcoulomb
    chern = σzτz[idx]
    ν = eachindex(ϵsorted) ./ length(ϵsorted)
    i = 1
    while (νF+4)/8 > ν[i]
        i += 1
    end
    if i<length(chern)
        # i = i-1
        ϵF = (ϵsorted[i+1] + ϵsorted[i])/2 
        i = length(ϵsorted[ϵsorted .<=hf.μ])
        Δ = ϵsorted[i+1] - ϵsorted[i]
    else
        ϵF = ϵsorted[end]
        Δ = 0 
    end
    return Δ
end

## plot Hartree Fock spectra collectively
function plot_spectra_collective(metadatas::Vector{String};savename::String="tmp.pdf",titlestr::String=" ")
    fig = figure(figsize=(2.4,2.8))
    ϵFs = Float64[]
    Δs = Float64[]
    for j in eachindex(metadatas) 
        metadata = metadatas[j]
        hf = load(metadata,"hf");
        νF = 8*round(Int,(hf.ν+4)/8*length(hf.ϵk)) / length(hf.ϵk)-4
        ϵk = hf.ϵk
        σzτz = hf.σzτz
        idx = sortperm(ϵk[:])
        ϵsorted = ϵk[idx] 
        chern = σzτz[idx]
        pl=scatter(ones(length(ϵsorted))*hf.p/hf.q,ϵsorted,c=chern,cmap="coolwarm",s=2,vmin=-1,vmax=1,marker=".")
        # plot(ones(length(ϵsorted[ϵsorted.<=hf.μ]))*hf.p/hf.q,ϵsorted[ϵsorted.<=hf.μ],"bo",markeredgecolor="none",markersize=2)
        # plot(ones(length(ϵsorted[ϵsorted.>hf.μ]))*hf.p/hf.q,ϵsorted[ϵsorted.>hf.μ],"o",c="gray",markeredgecolor="none",markersize=2)
        if j == length(metadatas)
            # colorbar(pl,shrink=0.8)
        end
        plot([hf.p/hf.q-0.01,hf.p/hf.q+0.01],[hf.μ,hf.μ],":",c="k",lw=0.5)
        
        ν = eachindex(ϵsorted) ./ length(ϵsorted)
        i = 1
        while (νF+4)/8 > ν[i]
            i += 1
        end
        # push!(ϵFs,(ϵsorted[i+1] + ϵsorted[i])/2) 
        # push!(Δs,max((ϵsorted[i+1] - ϵsorted[i]),(ϵsorted[i] - ϵsorted[i-1])))
        ν = eachindex(ϵsorted) ./ length(ϵsorted)
        if i<length(chern)
            # i = i-1
            ϵF = (ϵsorted[i+1] + ϵsorted[i])/2 
            i = length(ϵsorted[ϵsorted .<=hf.μ])
            Δ = ϵsorted[i+1] - ϵsorted[i]
        else
            ϵF = ϵsorted[end]
            Δ = 0 
        end
        push!(Δs,Δ)
        # push!(Δs,(ϵsorted[length(ϵsorted[ϵsorted .<=hf.μ]) +1] - ϵsorted[length(ϵsorted[ϵsorted .<=hf.μ])]))
        # axhline((ϵsorted[i+1] + ϵsorted[i])/2,ls=":",c="gray")
    end 
    # title(titlestr)
    # xticks(collect(0.1:0.2:0.5))
    xlim([0,0.55])
    ylim([-40,40])
    ylabel("E (meV)")
    xlabel(L"ϕ/ϕ_0")
    # ticklist = [1/2,1/3,2/7,1/4,1/5,1/6,1/8,1/10,1/14]
    # ticklistLabels= [L"$\frac{1}{2}$",L"$\frac{1}{3}$",
    #             L"$\frac{2}{7}$",L"$\frac{1}{4}$",L"$\frac{1}{5}$",
    #             L"$\frac{1}{6}$",L"$\frac{1}{8}$",L"$\frac{1}{10}$",L"$\frac{1}{14}$"]
    # ticklist = [1/4,1/8]
    # ticklistLabels= [L"$\frac{1}{4}$",L"$\frac{1}{8}$"]
    # xticks(ticklist,ticklistLabels)
    # title(titlestr)
    tight_layout()
    savefig(savename,dpi=600,transparent=true)
    display(fig)
    close(fig)

    # println("Gap sizes: ", Δs)
    return Δs
end


## plot error and energies under Hartree Fock iterations 
function plot_hf_iterations(fname::String)
    iter_err = load(fname,"iter_err")
    iter_energy = load(fname,"iter_energy")
    fig,ax = subplots(figsize=(5,3))
    ax.plot(iter_err,"r-",label="iter_err")
    # ax.set_ylim([0,2])
    ax.set_yscale("log")
    ax1 = ax.twinx()
    ax1.plot(iter_energy,"b-",label="iter_energy")
    ax.legend(loc="upper right")
    ax1.legend(loc="lower right")
    xlabel("Iter")
    ax.set_ylabel("Iter err")
    ax1.set_ylabel("Iter energy")
    tight_layout()
    # savefig(joinpath(fpath,"figures/tmp.pdf"),transparent=true)
    display(fig)
    close(fig)
    return nothing 
end

### density matrix analysis 
function plot_density_matrix_bm(fname::String;ik::Int=3)
    # plot at a given k point, 2qx4 x 2qx4 matrix 
    hf = load(fname,"hf");
    fig = figure(figsize=(6,6))
    P0 = view(hf.P,:,:,ik) + 0.5I
    pl = imshow(abs.(P0),vmin=0,vmax=1,origin="lower")
    # colorbar(pl)
    tight_layout()
    savefig("test.pdf")
    display(fig)
    close(fig)
    return nothing
end


### density matrix analysis 
function plot_density_matrix_bm_valley_spinv0(fname::String;ik::Int=1,savename::String="test.png")
    # plot at a given k point, 2qx4 x 2qx4 matrix 
    hf = load(fname,"hf");
    P0 = reshape(view(hf.P,:,:,ik)+0.5I,2hf.q,4,2hf.q,4);
    # P0 = reshape(reshape(sum(hf.P,dims=3),8hf.q,8hf.q)/size(hf.P,3)+0.5I,2hf.q,4,2hf.q,4);
    fig,ax = subplots(1,figsize=(2.3,2))
    # states = ["K↑","K'↑","K↓","K'↓"]
    pl = 0
    for r in 1:1, c in 1:1
        pl=ax.imshow(abs.(P0[:,r+2(c-1),:,r+2(c-1)]),extent=(1,2hf.q+1,1,2hf.q+1).-0.5,vmin=0,vmax=1,origin="lower",cmap="Reds")
        cbar = colorbar(pl,ax=ax,fraction=0.04, pad=0.1)
        cbar.set_ticks(collect(0:0.2:1))
        # colorbar(pl,ax=ax)
        # ax[r,c].set_title(states[r+2(c-1)])
        # ax[r,c].text(hf.q*0.94,2hf.q*0.9,states[r+2(c-1)],fontsize=12,color="k")
        # ax.set_xticks([1,hf.q,2hf.q])
        # ax.set_yticks([1,hf.q,2hf.q])
        # ax.set_xticklabels(["1","q","2q"])
        # ax.set_yticklabels(["1","q","2q"])
        ax.set_xticks([])
        ax.set_yticks([])
        # ax[r,c].axis("equal")
    end
    tight_layout()
    # subplots_adjust(hspace=0.02,wspace=0.05)
    # cbar_ax = fig.add_axes([0.74, 0.15, 0.15, 0.7])
    # cbar_ax.axis(false)
    # fig.colorbar(pl, ax=cbar_ax)
    
    savefig(savename,dpi=600,transparent=true)
    display(fig)
    close(fig)
    return nothing
end



### density matrix analysis 
function plot_density_matrix_bm_valley_spin(fname::String;ik::Int=1,savename::String="test.png")
    # plot at a given k point, 2qx4 x 2qx4 matrix 
    hf = load(fname,"hf");
    P0 = reshape(view(hf.P,:,:,ik)+0.5I,2hf.q,4,2hf.q,4);
    # P0 = reshape(reshape(sum(hf.P,dims=3),8hf.q,8hf.q)/size(hf.P,3)+0.5I,2hf.q,4,2hf.q,4);
    fig,ax = subplots(2,2,figsize=(3.2,3.2))
    states = ["K↑","K'↑","K↓","K'↓"]
    pl = 0
    for r in 1:2, c in 1:2 
        pl=ax[r,c].imshow(abs.(P0[:,r+2(c-1),:,r+2(c-1)]),extent=(1,2hf.q+1,1,2hf.q+1).-0.5,vmin=0,vmax=1,origin="lower",cmap="Blues")
        # colorbar(pl,ax=ax[r,c],fraction=0.046, pad=0.04)
        # ax[r,c].set_title(states[r+2(c-1)])
        ax[r,c].text(hf.q*0.94,2hf.q*0.9,states[r+2(c-1)],fontsize=12,color="k")
        ax[r,c].set_xticks([1,hf.q,2hf.q])
        ax[r,c].set_yticks([1,hf.q,2hf.q])
        ax[r,c].set_xticklabels(["1","q","2q"])
        ax[r,c].set_yticklabels(["1","q","2q"])
        # ax[r,c].axis("equal")
    end
    
    tight_layout()
    # subplots_adjust(hspace=0.02,wspace=0.05)
    # cbar_ax = fig.add_axes([0.74, 0.15, 0.15, 0.7])
    # cbar_ax.axis(false)
    # fig.colorbar(pl, ax=cbar_ax)
    
    savefig(savename,dpi=600,transparent=false)
    display(fig)
    close(fig)
    return nothing
end


### density matrix analysis 
function plot_density_matrix_bm_valley_spinv2(fname::String;ik::Int=1,savename::String="test.png")
    # plot at a given k point, 2qx4 x 2qx4 matrix 
    hf = load(fname,"hf");
    P0 = zeros(Float64,2hf.q,4,2hf.q,4)
    # for ik in 1:size(hf.P,3)
    for i in ik:ik
        P0 .+= abs.(reshape(view(hf.P,:,:,i)+0.5I,2hf.q,4,2hf.q,4));
    end
    P0 ./= 1 #size(hf.P,3)
    # P0 = reshape(reshape(sum(hf.P,dims=3),8hf.q,8hf.q)/size(hf.P,3)+0.5I,2hf.q,4,2hf.q,4);
    fig,ax = subplots(1,4,figsize=(8,2),sharex=true,sharey=true)
    states = ["K↑","K'↑","K↓","K'↓"]
    pl = 0
    for i in 1:4
        pl=ax[i].imshow(abs.(P0[:,i,:,i]),extent=(1,2hf.q+1,1,2hf.q+1).-0.5,vmin=0,vmax=1,origin="lower",cmap="coolwarm")
        if i==4
            colorbar(pl,ax=ax[i],fraction=0.046, pad=0.04)
        end
        # ax[r,c].set_title(states[r+2(c-1)])
        ax[i].text(hf.q*0.94,2hf.q*0.9,states[i],fontsize=14,color="w")
    end
    ax[1].set_xticks([1,hf.q,2hf.q])
    ax[1].set_yticks([1,hf.q,2hf.q])
    ax[1].set_xticklabels(["1","q","2q"])
    ax[1].set_yticklabels(["1","q","2q"])
    tight_layout()
    subplots_adjust(wspace=0.01)
    # cbar_ax = fig.add_axes([0.85, 0.15, 0.15, 0.7])
    # cbar_ax.axis(false)
    # fig.colorbar(pl, ax=cbar_ax)
    
    savefig(savename,dpi=500,transparent=true)
    display(fig)
    close(fig)
    return nothing
end


### density matrix analysis 
function plot_density_matrix_bm_valley_spinv3(fname::String;ik::Int=1,savename::String="test.png")
    # plot at a given k point, 2qx4 x 2qx4 matrix 
    hf = load(fname,"hf");
    # P0 = reshape(view(hf.P,:,:,ik)+0.5I,2hf.q,4,2hf.q,4);
    P0 = reshape(reshape(sum(hf.P,dims=3),8hf.q,8hf.q)/size(hf.P,3)+0.5I,2hf.q,4,2hf.q,4);
    fig = figure(figsize=(5,4))
    ax = fig.add_subplot(projection="3d")
    states = ["K↑","K'↑","K↓","K'↓"]
    cmaps = ["Blues","Reds","Greens","Purples"]
    pl = 0
    for i in 1:4
        coords = collect(1:2hf.q)
        Y = reshape(coords,:,1) .+ reshape(coords,1,:)*0.0
        X = reshape(coords,:,1)*0.0 .+ reshape(coords,1,:)
        # pl=ax.contourf(collect(1:2hf.q),collect(1:2hf.q),abs.(P0[:,i,:,i]),vmin=0,vmax=1,cmap="hot",zdir="z", offset=i-1)
        ax.plot_surface(X,Y, (i-1) * ones(size(P0,1),size(P0,3)), rstride=1, cstride=1, vmin=0,vmax=1,
                    facecolors=get_cmap("coolwarm")(abs.(P0[:,i,:,i])), shade=false, alpha=0.8)
        # colorbar(pl,ax=ax[r,c],fraction=0.046, pad=0.04)
        # ax[r,c].set_title(states[r+2(c-1)])
        # ax.text(hf.q*0.94,2hf.q*0.9,states[i],fontsize=14,color="w")
    end
    ax.set_xticks([1,hf.q,2hf.q])
    ax.set_yticks([1,hf.q,2hf.q])
    ax.set_zlim((-0.5,3.5))
    ax.set_xticklabels(["1","q","2q"])
    ax.set_yticklabels(["1","q","2q"])
    ax.set_zticks([0,1,2,3])
    ax.set_zticklabels(states)
    ax.axis(false)
    tight_layout()
    # subplots_adjust(wspace=0.01)
    # cbar_ax = fig.add_axes([0.85, 0.15, 0.15, 0.7])
    # cbar_ax.axis(false)
    # fig.colorbar(pl, ax=cbar_ax)
    
    savefig(savename,dpi=500,transparent=false)
    display(fig)
    close(fig)
    return nothing
end

function plot_density_matrix_global_order_parameters(fname::String)
    hf = load(fname,"hf");
    
    P = reshape(hf.P,8hf.q,8hf.q,hf.q,hf.nq,hf.nq)
    P = reshape(permutedims(P,(1,2,4,3,5)),8hf.q,8hf.q,:)
    PP = P .+ 0.5*reshape(Array{ComplexF64}(I,8hf.q,8hf.q),8hf.q,8hf.q,1)
    s0 = ComplexF64[1 0;0 1]
    sx = ComplexF64[0 1;1 0]
    sy = ComplexF64[0 -1im;1im 0]
    sz = ComplexF64[1 0;0 -1]
    Iq = Array{ComplexF64}(I,2hf.q,2hf.q)
    # Iq = diagm([(-1)^((hf.q-1)÷i) for i in 1:(2hf.q)])
    Iq = diagm([(-1)^(i) for i in 1:(2hf.q)])
    Os = kron(s0,kron(sz,Iq))
    Sk = reshape( sum(P.*reshape(Os,8hf.q,8hf.q,1),dims=(1,2)),hf.nq*hf.q,hf.nq) ./(hf.q)
    # println(sum(Sk)/length(Sk))
    fig = figure(figsize=(8,2))
    avgval = sum(real(Sk)) / length(Sk)
    δval = maximum(abs.(real(Sk).-avgval)) 
    println(δval)
    pl=imshow(real(Sk)',origin="lower",extent=(0,1,0,1/hf.q),cmap="bwr")
    colorbar(pl)
    xlabel(L"k_1")
    ylabel(L"k_2")
    title(L"s_z")
    tight_layout()
    # savefig("tmp.png",dpi=500,transparent=true)
    display(fig)
    close(fig)

    return nothing
end

# plot_density_matrix_global_order_parameters(metadata)



function plot_density_matrix_valley_spin_density_tL2(fname::String)
    hf = load(fname,"hf");
    P = reshape(hf.P,8hf.q,8hf.q,hf.q,hf.nq,hf.nq)
    P = reshape(permutedims(P,(1,2,4,3,5)),8hf.q,8hf.q,:)
    P = P .+ 0.5*reshape(Array{ComplexF64}(I,8hf.q,8hf.q),8hf.q,8hf.q,1)
    P = reshape(P,2hf.q,4,2hf.q,4,:)
    Iq = reshape(Array{ComplexF64}(I,2hf.q,2hf.q),2hf.q,2hf.q,1)
    strs = ["K↑","K'↑","K↓","K'↓"]
    fig, ax = subplots(4,1,figsize=(6,4))
    for i in 1:4 
        nk = reshape( sum(P[:,i,:,i,:].*Iq,dims=(1,2)) ./(hf.q), hf.q*hf.nq,hf.nq)
        pl=ax[i].imshow(real(nk)',origin="lower",extent=(0,1,0,1/hf.q),cmap="Reds",vmin=0,vmax=1.0)
        # colorbar(pl,ax=ax[i])
    end
    for i in 1:4
        ax[4].set_xlabel(L"k_1")
        if i!=4 
            ax[i].set_xticklabels([])
        end
        ax[i].set_ylabel(strs[i])
    end
    tight_layout()
    subplots_adjust(hspace=-0.1,wspace=0)
    savefig("tmp_$(hf.p)_$(hf.q).png",dpi=500,transparent=true)
    display(fig)
    close(fig)

    return nothing
end

# plot_density_matrix_valley_spin_density_tL2(metadata)

### strong coupling basis 
function plot_density_matrix_strong_coupling(fname::String,fname0::String)
    ik = 1
    hf = load(fname,"hf");
    hf0 = load(fname0,"hf");
    H0 = hf0.H
    P = hf.P 
    tmpH0 = view(H0,:,:,ik)
    P0 = view(P,:,:,ik)+0.5I
    F = eigen(Hermitian(tmpH0))
    Pstrong = transpose(F.vectors) * P0 * conj.(F.vectors )

    fig = figure(figsize=(6,6))
    pl=imshow(abs.(Pstrong),vmin=0,vmax=1,origin="lower")
    colorbar(pl,fraction=0.046, pad=0.04)
    tight_layout()
    display(fig)
    close(fig)
    return nothing 
end


### strong coupling basis  valley spin
function plot_density_matrix_strong_coupling_valley_spin(fname::String,fname0::String)
    ik = 1
    hf = load(fname,"hf");
    hf0 = load(fname0,"hf");
    H0 = hf0.H
    P = hf.P 
    tmpH0 = reshape(view(H0,:,:,ik),2hf.q,4,2hf.q,4)
    P0 = reshape(view(P,:,:,ik)+0.5I,2hf.q,4,2hf.q,4)
    Pstrong = zeros(ComplexF64,2hf.q,2hf.q,4)
    for iηs in 1:4
        F = eigen(Hermitian(tmpH0[:,iηs,:,iηs]))
        vec = F.vectors
        Pstrong[:,:,iηs] = transpose(vec) * P0[:,iηs,:,iηs] * conj.(vec) 
    end 
    fig,ax = subplots(2,2,figsize=(6,6))
    states = ["K↑","K'↑","K↓","K'↓"]
    for r in 1:2, c in 1:2 
        pl=ax[r,c].imshow(abs.(Pstrong[:,:,r+2(c-1)]),vmin=0,vmax=1,origin="lower")
        colorbar(pl,ax=ax[r,c],fraction=0.046, pad=0.04)
        ax[r,c].set_title(states[r+2(c-1)])
    end
    tight_layout()
    savefig("test.pdf")
    display(fig)
    close(fig)
    return nothing 
end


### strong coupling basis  valley spin
function plot_density_matrix_sublattice(fname::String)
    ik = 1
    hf = load(fname,"hf");
    H0 = hf.Σz0
    P = hf.P 
    tmpH0 = reshape(view(H0,:,:,ik),2hf.q,4,2hf.q,4)
    P0 = reshape(view(P,:,:,ik)+0.5I,2hf.q,4,2hf.q,4)
    Pstrong = zeros(ComplexF64,2hf.q,2hf.q,4)
    for iηs in 1:4
        F = eigen(Hermitian(tmpH0[:,iηs,:,iηs]))
        vec = F.vectors
        Pstrong[:,:,iηs] = transpose(vec) * P0[:,iηs,:,iηs] * conj.(vec) 
    end 
    fig,ax = subplots(2,2,figsize=(6,6))
    states = ["K↑","K'↑","K↓","K'↓"]
    for r in 1:2, c in 1:2 
        pl=ax[r,c].imshow(abs.(Pstrong[:,:,r+2(c-1)]),vmin=0,vmax=1,origin="lower")
        colorbar(pl,ax=ax[r,c],fraction=0.046, pad=0.04)
        ax[r,c].set_title(states[r+2(c-1)])
    end
    tight_layout()
    savefig("test.pdf")
    display(fig)
    close(fig)
    return nothing 
end




### strong coupling basis  valley spin
function plot_density_matrix_sublattice_full(fname::String)
    ik = 1
    hf = load(fname,"hf");
    H0 = hf.Σz0
    P = hf.P 
    tmpH0 = reshape(view(H0,:,:,ik),2hf.q,4,2hf.q,4)
    P0 = reshape(view(P,:,:,ik)+0.5I,2hf.q,4,2hf.q,4)
    Pstrong = zeros(ComplexF64,2hf.q,4,2hf.q,4)
    for iηs in 1:4, iηs1 in 1:4
        F = eigen(Hermitian(tmpH0[:,iηs,:,iηs]))
        vec = F.vectors
        F1 = eigen(Hermitian(tmpH0[:,iηs1,:,iηs1]))
        vec1 = F1.vectors
        Pstrong[:,iηs,:,iηs1] = transpose(vec) * P0[:,iηs,:,iηs1] * conj.(vec1) 
    end 
    fig = figure(figsize=(6,6))
    pl=imshow(abs.(reshape(Pstrong,8hf.q,8hf.q)),vmin=0,vmax=1,origin="lower")
    colorbar(pl,fraction=0.046, pad=0.04)
    tight_layout()
    savefig("test.pdf")
    display(fig)
    close(fig)
    return nothing 
end