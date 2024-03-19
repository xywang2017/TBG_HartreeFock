using PyPlot 

## plot Hartree Fock spectra
function plot_spectra(metadata::String;savename::String="tmp.pdf")
    hf = load(metadata,"hf");
    ϵk = hf.ϵk
    σzτz = hf.σzτz

    fig = figure(figsize=(3,3))
    idx = sortperm(ϵk[:])
    ϵsorted = ϵk[idx] #./Vcoulomb
    chern = σzτz[idx]

    pl=scatter(ones(length(ϵsorted))*hf.p/hf.q,ϵsorted,c=chern,cmap="coolwarm",s=8,vmin=-1,vmax=1)
    colorbar(pl)
    
    axhline(hf.μ,ls=":",c="gray")
    ylabel("E (meV)")
    xlabel(L"ϕ/ϕ_0")
    tight_layout()
    savefig(savename,transparent=true,dpi=600)
    display(fig)
    close(fig)
    return nothing
end


## plot Hartree Fock spectra
function plot_spectrav2(metadata::String;savename::String="tmp.pdf")
    hf = load(metadata,"hf");
    ϵk = zeros(Float64,size(hf.H,1),size(hf.H,3))
    H = reshape(hf.H,2hf.q,4,2hf.q,4,:) 
    H1 = reshape(H[:,[1,4],:,[1,4],:],4hf.q,4hf.q,:)
    H2 = reshape(H[:,[2,3],:,[2,3],:],4hf.q,4hf.q,:)


    a0 = Float64[1 0; 0 1]
    az = Float64[1 0; 0 -1]
    ηz = zeros(Float64,size(ϵk))
    Oηz = kron(a0,kron(az,Array{Float64,2}(I,2hf.q,2hf.q)))
    tmpOηz = reshape(Oηz,2hf.q,4,2hf.q,4)
    O1 = reshape(tmpOηz[:,[1,4],:,[1,4]],4hf.q,4hf.q)
    O2 = reshape(tmpOηz[:,[2,3],:,[2,3]],4hf.q,4hf.q)

    for ik in 1:size(ϵk,2)
        F = eigen(Hermitian(H1[:,:,ik]))
        ϵk[1:(4hf.q),ik] = real(F.values )
        for iq in 1:size(H1,1)
            ηz[iq,ik] = real(F.vectors[:,iq]'*O1*F.vectors[:,iq] )
        end
        F = eigen(Hermitian(H2[:,:,ik]))
        ϵk[(4hf.q+1):end,ik] = real(F.values )
        for iq in 1:size(H2,1)
            ηz[iq+4hf.q,ik] = real(F.vectors[:,iq]'*O2*F.vectors[:,iq] )
        end
    end
    
    if !occursin("_tL_",metadata)
        ϵk = reshape(ϵk,8hf.q,hf.q,hf.nq^2)[:,1,:]
        ηz = reshape(ηz,8hf.q,hf.q,hf.nq^2)[:,1,:]
    end
    
    fig = figure(figsize=(3,3))
    pl=scatter(ones(length(ϵk[ϵk .< hf.μ]))*hf.p/hf.q,ϵk[ϵk .< hf.μ],c=abs.(ηz[ϵk .< hf.μ]),cmap="coolwarm",s=4,vmin=0,vmax=1,marker="o")
    # scatter(ones(length(ϵk[ϵk .>= hf.μ]))*hf.p/hf.q,ϵk[ϵk .>= hf.μ],c="gray",s=2,marker="o")
    pl=scatter(ones(length(ϵk[ϵk .> hf.μ]))*hf.p/hf.q,ϵk[ϵk .> hf.μ],c=abs.(ηz[ϵk .> hf.μ]),cmap="coolwarm",s=4,vmin=0,vmax=1,marker="o")
    colorbar(pl)
    axhline(hf.μ,ls=":",c="gray")
    ylabel("E (meV)")
    xlabel(L"ϕ/ϕ_0")
    tight_layout()
    savefig(savename,transparent=true,dpi=600)
    display(fig)
    close(fig)
    return nothing
end



## plot Hartree Fock spectra
function plot_spectrav3(metadata::String;savename::String="tmp.pdf")
    hf = load(metadata,"hf");
    H = reshape(hf.H,2hf.q,4,2hf.q,4,:,hf.nq^2)
    H = H[:,:,:,:,1,:] 
    Σz0 = reshape(hf.Σz0,2hf.q,4,2hf.q,4,:)
    subblocks = [[1],[2],[3,4]]
    lbls = ["K↑","K'↑","{K↓,K'↓}"]
    ϵk = []
    σz = []
    for sb in subblocks
        ϵk0 = zeros(Float64,2hf.q*length(sb),hf.nq^2)
        σz0 = zeros(Float64,2hf.q*length(sb),hf.nq^2)
        for ik in 1:size(ϵk0,2)
            F = eigen(Hermitian(reshape(H[:,sb,:,sb,ik],:,length(sb)*2hf.q)))
            ϵk0[:,ik] = F.values 
            for iq in 1:size(ϵk0,1)
                σz0[iq,ik] = real(F.vectors[:,iq]'*reshape(Σz0[:,sb,:,sb,ik],:,length(sb)*2hf.q)*F.vectors[:,iq])
            end
        end
        push!(ϵk,ϵk0[:])
        push!(σz,σz0[:])
    end

    fig = figure(figsize=(3,3))
    for ib in 1:length(subblocks)
        ϵ0 = ϵk[ib]
        σ0 = σz[ib]
        pl = scatter(ones(length(ϵ0))*ib,ϵ0,c=σ0,cmap="coolwarm",s=6,vmin=-0.4,vmax=0.4,marker="o")
        if ib ==1 
            colorbar(pl)
        end
    end
    axhline(hf.μ,ls=":",c="gray")
    ylabel("E (meV)")
    xlim([0.5,3.5])
    xticks(collect(eachindex(subblocks)),lbls)
    tight_layout()
    savefig(savename,transparent=true,dpi=600)
    display(fig)
    close(fig)
    return nothing
end

## plot Hartree Fock spectra per flavor
function plot_spectra_flavor(metadata::String;savename::String="tmp.pdf")
    hf = load(metadata,"hf");
    νF = 8*round(Int,(hf.ν+4)/8*length(hf.ϵk)) / length(hf.ϵk)-4
    ϵk = hf.ϵk
    σzτz = hf.σzτz
    params = hf.params 

    H = reshape(hf.H,2hf.q,4,2hf.q,4,:)
    Σz = reshape(hf.Σz0,2hf.q,4,2hf.q,4,:)
    ϵk = zeros(Float64,2hf.q)
    σz = zeros(Float64,2hf.q)
    colors=["tab:blue","tab:red","tab:green","tab:purple"]
    str = ["K↑","K'↑","K↓","K'↓"]

    fig, ax = subplots(2,2,figsize=(4,4),sharex=true,sharey=true)
    for r in 1:2, c in 1:2 
        elem = (r-1)*2 + c
        for ik in 1:size(H,5) 
            tmpH = view(H,:,elem,:,elem,ik)
            F = eigen(Hermitian(tmpH))
            ϵk .= real(F.values)
            for j in eachindex(σz)
                σz[j] = real(F.vectors[:,j]'*view(Σz,:,elem,:,elem,ik)*F.vectors[:,j])
            end
            # ax[r,c].plot(ones(length(ϵk))*hf.p/hf.q,ϵk,"o",c=colors[elem],ms=2)
            ax[r,c].scatter(ones(length(ϵk))*hf.p/hf.q,ϵk,c=σz,cmap="coolwarm",s=4,vmin=-0.5,vmax=0.5)
        end
        ax[r,c].axhline(hf.μ,ls=":",c="gray")
        if r==2 
            ax[r,c].set_xlabel(L"ϕ/ϕ_0")
        end
        if c==1 
            ax[r,c].set_ylabel("E (meV)")
        end
        ax[r,c].text(0.1,20,str[elem])
    end
    tight_layout()
    savefig("1.20_degrees_-3_-1_1_5.png",dpi=500)
    display(fig)
    close(fig)
    return nothing
end

function plot_spectra_flavorv1(metadata::String;savename::String="tmp.pdf")
    hf = load(metadata,"hf");
    νF = 8*round(Int,(hf.ν+4)/8*length(hf.ϵk)) / length(hf.ϵk)-4
    ϵk = hf.ϵk
    σzτz = hf.σzτz
    params = hf.params 

    H = reshape(hf.H,2hf.q,4,2hf.q,4,:)
    Σz = reshape(hf.Σz0,2hf.q,4,2hf.q,4,:)
    ϵk = zeros(Float64,2hf.q)
    σz = zeros(Float64,2hf.q)
    colors=["tab:blue","tab:red","tab:green","tab:purple"]
    str = ["K↑","K'↑","K↓","K'↓"]

    fig, ax = subplots(1,1,figsize=(3,2.5))
    pl = 0
    for elem in 1:4
        for ik in 1:size(H,5) 
            tmpH = view(H,:,elem,:,elem,ik)
            F = eigen(Hermitian(tmpH))
            ϵk .= real(F.values)
            for j in eachindex(σz)
                σz[j] = real(F.vectors[:,j]'*view(Σz,:,elem,:,elem,ik)*F.vectors[:,j])
            end
            # ax[r,c].plot(ones(length(ϵk))*hf.p/hf.q,ϵk,"o",c=colors[elem],ms=2)
            pl =ax.scatter(ones(length(ϵk))*elem,ϵk,c=σz,cmap="coolwarm",s=6,vmin=-1,vmax=1)
        end 
    end
    # colorbar(pl,ax=ax)
    ax.axhline(hf.μ,ls=":",c="gray")
    ax.set_xlim([0,5])
    ax.set_ylim(-35,39)
    ax.set_ylabel("E (meV)")
    ax.set_xticks([1,2,3,4])
    ax.set_xticklabels(str)
    tight_layout()
    savefig("tmp.png",dpi=500,transparent=true)
    display(fig)
    close(fig)
    return nothing
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
    Pz = sum(chern[1:i])/ i
    return Δ,Pz 
end

## plot Hartree Fock spectra collectively
function plot_spectra_collective(metadatas::Vector{String};savename::String="tmp.pdf",titlestr::String=" ",indices::Vector{Int}=Int[])
    fig = figure(figsize=(3,2.5))
    # fig = figure(figsize=(5,4))
    ϕs = Float64[]
    Δs = Float64[]
    cmap =["coolwarm","bwr"]
    μs = Float64[]
    for j in eachindex(metadatas) 
        metadata = metadatas[j]
        hf = load(metadata,"hf");
        νF = 8*round(Int,(hf.ν+4)/8*length(hf.ϵk)) / length(hf.ϵk)-4
        # if true # (hf.p/hf.q!=1/3 && hf.p/hf.q!=4/11 && hf.p/hf.q!=1/2)
        # if hf.p/hf.q < 1/3
        ϵk = hf.ϵk
        σzτz = hf.σzτz
        idx = sortperm(ϵk[:])
        ϵsorted = ϵk[idx] 
        chern = σzτz[idx]
        pl=scatter(ones(length(ϵsorted))*hf.p/hf.q,ϵsorted,c=chern,cmap="coolwarm",s=3,vmin=-1,vmax=1,marker="o",edgecolor="none")
        # if j in indices
        #     plot(ones(length(ϵsorted[ϵsorted.<=hf.μ]))*hf.p/hf.q,ϵsorted[ϵsorted.<=hf.μ],"o",c=[0,0.6,0],markeredgecolor="none",markersize=1.5)
        # else
        #     plot(ones(length(ϵsorted[ϵsorted.<=hf.μ]))*hf.p/hf.q,ϵsorted[ϵsorted.<=hf.μ],"o",c=[0.8,0,0],markeredgecolor="none",markersize=1.5)
        # end
        # plot(ones(length(ϵsorted[ϵsorted.>hf.μ]))*hf.p/hf.q,ϵsorted[ϵsorted.>hf.μ],"o",c="gray",markeredgecolor="none",markersize=1.5)
        if j == length(metadatas)
             colorbar(pl,shrink=0.8)
        end
        # plot([hf.p/hf.q-0.01,hf.p/hf.q+0.01],[hf.μ,hf.μ],":",c="k",lw=0.5)
        push!(μs,hf.μ)
        push!(ϕs,hf.p/hf.q)
        
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
        # end
        push!(Δs,Δ)
    end
        # push!(Δs,(ϵsorted[length(ϵsorted[ϵsorted .<=hf.μ]) +1] - ϵsorted[length(ϵsorted[ϵsorted .<=hf.μ])]))
        # axhline((ϵsorted[i+1] + ϵsorted[i])/2,ls=":",c="gray")
    end 
    plot(ϕs,μs,":",c="k",lw=0.5)
    # title(titlestr)
    # xticks(collect(0.1:0.2:0.5))
    # axvline([(2//9+1//4)/2],ls=":",c="k")
    # axvline([1//8+1//7]/2,ls=":",c="k")
    xlim([0,0.55])
    xticks([0.1,0.3,0.5])
    ylim([-44,44])
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



## plot Hartree Fock spectra collectively, color by spin + valley polarization
function plot_spectra_collectivev2(metadatas::Vector{String};savename::String="tmp.pdf",titlestr::String=" ",indices::Vector{Int}=Int[])
    # fig = figure(figsize=(3,2.5))
    fig = figure(figsize=(3.3,3))
    # fig = figure(figsize=(6,4))
    ϕs = Float64[]
    μs = Float64[]
    a0 = Float64[1 0; 0 1]
    az = Float64[1 0; 0 -1]
    for j in eachindex(metadatas) 
        metadata = metadatas[j]
        hf = load(metadata,"hf");
        ϵk = zeros(Float64,size(hf.H,1),size(hf.H,3))
        sz = zeros(Float64,size(ϵk))
        ηz = zeros(Float64,size(ϵk))
        Osz = kron(az,kron(a0,Array{Float64,2}(I,2hf.q,2hf.q)))
        Oηz = kron(a0,kron(az,Array{Float64,2}(I,2hf.q,2hf.q)))
        for ik in 1:size(ϵk,2)
            H = hf.H[:,:,ik]
            H[abs.(H) .<1e-2] .= 0.0
            F = eigen(Hermitian(H))
            ϵk[:,ik] = F.values 
            for iq in 1:size(hf.H,1)
                sz[iq,ik] = real(F.vectors[:,iq]'*Osz*F.vectors[:,iq])
                ηz[iq,ik] = real(F.vectors[:,iq]'*Oηz*F.vectors[:,iq] )
            end
        end
        if !occursin("_tL_",metadata)
            ϵk = reshape(ϵk,8hf.q,hf.q,hf.nq^2)[:,1,:]
            ηz = reshape(ηz,8hf.q,hf.q,hf.nq^2)[:,1,:]
        end
        pl=scatter(ones(length(ϵk[ϵk .< hf.μ]))*hf.p/hf.q,ϵk[ϵk .< hf.μ],c=abs.(ηz[ϵk .< hf.μ]),cmap="coolwarm",s=6,vmin=0,vmax=1,marker="o",edgecolors="none")
        scatter(ones(length(ϵk[ϵk .>= hf.μ]))*hf.p/hf.q,ϵk[ϵk .>= hf.μ],c="gray",s=6,marker="o",edgecolors="none")
        # pl=scatter(ones(length(ϵk[ϵk .> hf.μ]))*hf.p/hf.q,ϵk[ϵk .> hf.μ],c=abs.(ηz[ϵk .> hf.μ]),cmap="coolwarm",s=2,vmin=0,vmax=1,marker="o")
        if j == length(metadatas)
            #  colorbar(pl,shrink=0.8)
        end
        push!(μs,hf.μ)
        push!(ϕs,hf.p/hf.q)

        # if hf.q == 12  && hf.p == 1
        #     println(abs.(ηz[ϵk .< hf.μ]))
        # end
    end 
    plot(ϕs,μs,":",c="k",lw=0.5)
    xlim([0.01,0.55])
    ylim([-44,44])
    yticks(collect(-40:20:20),fontsize=13)
    xticks([0.1,0.2,0.3,0.4,0.5],fontsize=13)
    ylabel("E (meV)",fontsize=13)
    xlabel(L"ϕ/ϕ_0",fontsize=13)
    tight_layout()
    savefig(savename,dpi=600,transparent=true)
    display(fig)
    close(fig)
    return nothing
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
function plot_density_matrix_bm_half(fname::String;ik::Int=1,savename::String="test.png")
    # plot at a given k point, 2qx4 x 2qx4 matrix 
    hf = load(fname,"hf");
    fig = figure(figsize=(3.5,3.5))
    P0 = view(hf.P,(4hf.q+1):(8hf.q),(4hf.q+1):(8hf.q),ik) + 0.5I
    pl = imshow(abs.(P0),vmin=0,vmax=1,origin="lower",cmap="Blues",extent=(1,4hf.q+1,1,4hf.q+1).-0.5)
    xticks([])
    yticks([])
    axhline(2hf.q+0.5,ls=":",c="gray")
    axvline(2hf.q+0.5,ls=":",c="gray")
    colorbar(pl,shrink=0.8)
    tight_layout()
    savefig(savename,dpi=600,transparent=true)
    display(fig)
    close(fig)
    return nothing
end

function plot_density_matrix_bm(fname::String;ik::Int=1,savename::String="test.png")
    # plot at a given k point, 2qx4 x 2qx4 matrix 
    hf = load(fname,"hf");
    # fig = figure(figsize=(3.6,3))
    fig = figure(figsize=(2.5,2.5))
    P0 = view(hf.P,:,:,ik) + 0.5I
    # P0 = view(hf.H,:,:,ik)
    pl = imshow(abs.(P0),vmin=0,vmax=1,origin="lower",cmap="Blues",extent=(1,8hf.q+1,1,8hf.q+1).-0.5)
    # pl = imshow(abs.(P0),vmin=0,origin="lower",cmap="Blues",extent=(1,8hf.q+1,1,8hf.q+1).-0.5)
    xticks([])
    yticks([])
    for r in [2hf.q,4hf.q,6hf.q]
        axhline(r+0.5,ls=":",c="gray")
        axvline(r+0.5,ls=":",c="gray")
    end
    colorbar(pl,shrink=0.9)
    tight_layout()
    savefig(savename,dpi=600,transparent=true)
    display(fig)
    close(fig)
    return nothing
end

### density matrix analysis 
function plot_density_matrix_bm_valley_spinv0(fname::String;ik::Int=1,savename::String="test.png",jj::Int=1,titlestr::String="")
    # plot at a given k point, 2qx4 x 2qx4 matrix 
    hf = load(fname,"hf");
    P0 = reshape(view(hf.P,:,:,ik)+0.5I,2hf.q,4,2hf.q,4);
    # P0 = reshape(reshape(sum(hf.P,dims=3),8hf.q,8hf.q)/size(hf.P,3)+0.5I,2hf.q,4,2hf.q,4);
    fig,ax = subplots(1,figsize=(2.5,2.5))
    # states = ["K↑","K'↑","K↓","K'↓"]
    pl = 0
    
    pl=ax.imshow(abs.(P0[:,jj,:,jj]),extent=(1,2hf.q+1,1,2hf.q+1).-0.5,vmin=0,vmax=1,origin="lower",cmap="Reds")
    # cbar = colorbar(pl,ax=ax,fraction=0.04, pad=0.1)
    # cbar.set_ticks(collect(0:0.2:1))
    # colorbar(pl,ax=ax)
    # ax[r,c].set_title(states[r+2(c-1)])
    # ax[r,c].text(hf.q*0.94,2hf.q*0.9,states[r+2(c-1)],fontsize=12,color="k")
    # ax.set_xticks([1,hf.q,2hf.q])
    # ax.set_yticks([1,hf.q,2hf.q])
    # ax.set_xticklabels(["1","q","2q"])
    # ax.set_yticklabels(["1","q","2q"])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(titlestr)
    # ax[r,c].axis("equal")
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
        # ax[r,c].text(hf.q*0.94,2hf.q*0.9,states[r+2(c-1)],fontsize=12,color="k")
        ax[r,c].set_xticks([1,hf.q,2hf.q])
        ax[r,c].set_yticks([1,hf.q,2hf.q])
        ax[r,c].set_xticklabels(["1","q","2q"])
        ax[r,c].set_yticklabels(["1","q","2q"])
        ax[r,c].axis("off")
    end
    
    tight_layout()
    subplots_adjust(hspace=0.02,wspace=0.05)
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
        # ax[i].text(hf.q*0.94,2hf.q*0.9,states[i],fontsize=14,color="w")
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
    cmaps = ["Blues","Blues","Blues","Reds"]
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
    
    P = reshape(hf.P,8hf.q,8hf.q,:,hf.nq,hf.nq)
    P = reshape(permutedims(P,(1,2,4,3,5)),8hf.q,8hf.q,:)
    PP = P .+ 0.5*reshape(Array{ComplexF64}(I,8hf.q,8hf.q),8hf.q,8hf.q,1)
    s0 = ComplexF64[1 0;0 1]
    sx = ComplexF64[0 1;1 0]
    sy = ComplexF64[0 -1im;1im 0]
    sz = ComplexF64[1 0;0 -1]
    Iq = Array{ComplexF64}(I,2hf.q,2hf.q)
    Iq = diagm([(-1)^((hf.q-1)÷i) for i in 1:(2hf.q)])
    # Iq = diagm([(-1)^(i) for i in 1:(2hf.q)])
    Os = kron(s0,kron(sz,Iq))
    Sk = reshape( [tr(transpose(P[:,:,ik])*Os)/hf.q for ik in 1:size(P,3)], (:,hf.nq) )
    Sk = reshape( P[7hf.q,7hf.q,:].+0.5, (:,hf.nq) ) 
    # println(sum(Sk)/length(Sk))
    fig = figure(figsize=(6,1.5))
    avgval = sum(real(Sk)) / length(Sk)
    δval = maximum(abs.(real(Sk).-avgval)) 
    println(δval)
    pl=imshow(real(Sk)',origin="lower",extent=(0,1,0,1/hf.q),cmap="coolwarm",vmin=0.9*minimum(abs.(real(Sk))),vmax=1.1*maximum(abs.(real(Sk))))
    # Coords = reshape(1:(hf.nq*hf.q),:,1)/(hf.nq*hf.q) .+ 1im*reshape(1:hf.nq,1,:)/(hf.nq*hf.q)
    # pl = scatter(real(Coords),imag(Coords),s=100,c=real(Sk),cmap="viridis",vmin=0,vmax=1.2*maximum(abs.(real(Sk))))
    colorbar(pl)
    # axhline(1/(hf.nq*hf.q),c="k",ls=":")
    for line in 1:hf.q 
        axvline(line/(hf.q),c="k",ls=":")
    end
    xticks([])
    yticks([])
    # xlabel(L"k_1")
    # ylabel(L"k_2")
    # title(L"s_z")
    # title("")
    tight_layout()
    savefig("tmp.png",dpi=500,transparent=true)
    display(fig)
    close(fig)

    return nothing
end

# plot_density_matrix_global_order_parameters(metadata)


function test_tL2_breaking(fname::String)
    hf = load(fname,"hf");
    P = reshape(hf.P,8hf.q,8hf.q,:,hf.nq^2)
    fluctP = P .- reshape(sum(P,dims=3),(8hf.q,8hf.q,1,:))/size(P,3)
    normfluctP = norm(fluctP) /length(fluctP)
    if normfluctP > 1e-4 
        println(normfluctP)
    else 
        println("tL2 preserved with resolution: ",normfluctP)
    end
    return nothing
end


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
        pl=ax[i].imshow(real(nk)',origin="lower",extent=(0,1,0,1/hf.q),cmap="Blues",vmin=0,vmax=1.0)
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
        pl=ax[r,c].imshow(abs.(Pstrong[:,:,r+2(c-1)]),vmin=0,vmax=1,origin="lower",cmap="hot")
        colorbar(pl,ax=ax[r,c],fraction=0.046, pad=0.04)
        ax[r,c].set_title(states[r+2(c-1)])
    end
    tight_layout()
    savefig("test.png",dpi=500,transparent=true)
    display(fig)
    close(fig)
    return nothing 
end


### strong coupling basis  valley spin
function plot_density_matrix_strong_coupling_valley_spin_v1(fname::String,fname0::String;ik::Int=1)
    hf = load(fname,"hf");
    hf0 = load(fname0,"hf");
    H0 = hf0.H
    P = hf.P 
    tmpH0 = reshape(view(H0,:,:,ik),2hf.q,4,2hf.q,4)
    P0 = reshape(view(P,:,:,ik)+0.5I,2hf.q,4,2hf.q,4)
    Pstrong = zeros(ComplexF64,2hf.q,4,2hf.q,4)
    vec = zeros(ComplexF64,2hf.q,2hf.q,4)
    for iηs in 1:4
        F = eigen(Hermitian(tmpH0[:,iηs,:,iηs]))
        vec[:,:,iηs] = F.vectors
    end 
    for iη1 in 1:4, iη2 in 1:4 
        Pstrong[:,iη1,:,iη2] = transpose(vec[:,:,iη1]) * P0[:,iη1,:,iη2] * conj.(vec[:,:,iη2]) 
    end
    fig = figure(figsize=(3.6,3))
    pl=imshow(abs.(reshape(Pstrong,8hf.q,:)),vmin=0,vmax=1,origin="lower",extent=(1,8hf.q+1,1,8hf.q+1).-0.5)
    
    xticks([])
    yticks([])
    for r in [2hf.q,4hf.q,6hf.q]
        axhline(r+0.5,ls=":",c="gray")
        axvline(r+0.5,ls=":",c="gray")
    end
    colorbar(pl,shrink=0.9)
    # fig = figure(figsize=(2.4,2))
    # plot(1:(hf.q+hf.p),σzτz0[1:hf.q+hf.p],"bo",ms=3)
    # plot((hf.q+hf.p+1):(2hf.q),σzτz0[(hf.q+hf.p+1):(2hf.q)],"ro",ms=3)
    # xticks([1,hf.q+hf.p,2hf.q],["1","q+p","2q"])
    # ylabel(L"\rm eigvals\ of\ σ_zτ_z")
    tight_layout()
    savefig("test.png",dpi=500,transparent=true)
    display(fig)
    close(fig)
    return nothing 
end




### strong coupling basis  valley spin
function plot_density_matrix_strong_coupling_valley_spin_v1_half(fname::String,fname0::String)
    ik = 1
    hf = load(fname,"hf");
    hf0 = load(fname0,"hf");
    H0 = hf0.H
    P = hf.P 
    tmpH0 = reshape(view(H0,:,:,ik),2hf.q,4,2hf.q,4)
    P0 = reshape(view(P,:,:,ik)+0.5I,2hf.q,4,2hf.q,4)
    Pstrong = zeros(ComplexF64,2hf.q,4,2hf.q,4)
    vec = zeros(ComplexF64,2hf.q,2hf.q,4)
    for iηs in 1:4
        F = eigen(Hermitian(tmpH0[:,iηs,:,iηs]))
        vec[:,:,iηs] = F.vectors
    end 
    for iη1 in 1:4, iη2 in 1:4 
        Pstrong[:,iη1,:,iη2] = transpose(vec[:,:,iη1]) * P0[:,iη1,:,iη2] * conj.(vec[:,:,iη2]) 
    end
    fig = figure(figsize=(2.5,2.5))
    pl=imshow(abs.(reshape(Pstrong,8hf.q,:))[(6hf.q+1):end,(6hf.q+1):end],vmin=0,vmax=1,origin="lower",cmap="hot")
    colorbar(pl,shrink=0.6)
    # axhline(2hf.q-0.5,c="w",ls=":")
    # axvline(2hf.q-0.5,c="w",ls=":")
    xticks([])
    yticks([])
    tight_layout()
    savefig("test.png",dpi=500,transparent=true)
    display(fig)
    close(fig)
    return nothing 
end


### strong coupling basis  valley spin
function plot_density_matrix_sublattice_v0(fname::String;titlestr::String="")
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
        # Pstrong[:,:,iηs] = vec' * P0[:,iηs,:,iηs] * vec
    end 
    fig,ax = subplots(figsize=(3,2.8))
    pl=ax.imshow(abs.(Pstrong[:,:,3]),vmin=0,vmax=1,origin="lower",cmap="hot")
    colorbar(pl,ax=ax,fraction=0.046, pad=0.04)
    ax.set_title(titlestr)
    tight_layout()
    savefig("figures/test_$(titlestr).png",dpi=500)
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
        # Pstrong[:,:,iηs] = vec' * P0[:,iηs,:,iηs] * vec
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
function plot_density_matrix_sublattice_full(fname::String;ik::Int=1)
    hf = load(fname,"hf");
    H0 = hf.Σz0
    P = hf.P 
    tmpH0 = reshape(view(H0,:,:,ik),2hf.q,4,2hf.q,4)
    σzτz0 = eigvals(Hermitian(tmpH0[:,1,:,1]))
    P0 = reshape(view(P,:,:,ik)+0.5I,2hf.q,4,2hf.q,4)
    Pstrong = zeros(ComplexF64,2hf.q,4,2hf.q,4)
    for iηs in 1:4, iηs1 in 1:4
        F = eigen(Hermitian(tmpH0[:,iηs,:,iηs]))
        vec = F.vectors
        F1 = eigen(Hermitian(tmpH0[:,iηs1,:,iηs1]))
        vec1 = F1.vectors
        Pstrong[:,iηs,:,iηs1] = transpose(vec) * P0[:,iηs,:,iηs1] * conj.(vec1) 
    end 
    fig = figure(figsize=(4.65,4))
    pl=imshow(abs.(reshape(Pstrong,8hf.q,8hf.q)),vmin=0,vmax=1,origin="lower",extent=(1,8hf.q+1,1,8hf.q+1).-0.5)
    xticks([])
    yticks([])
    for r in [2hf.q,4hf.q,6hf.q]
        axhline(r+0.5,ls=":",c="gray")
        axvline(r+0.5,ls=":",c="gray")
    end
    colorbar(pl,shrink=0.8)
    # fig = figure(figsize=(2.4,2))
    # plot(1:(hf.q+hf.p),σzτz0[1:hf.q+hf.p],"bo",ms=3)
    # plot((hf.q+hf.p+1):(2hf.q),σzτz0[(hf.q+hf.p+1):(2hf.q)],"ro",ms=3)
    # xticks([1,hf.q+hf.p,2hf.q],["1","q+p","2q"])
    # ylabel(L"\rm eigvals\ of\ σ_zτ_z")
    tight_layout()
    savefig("test.png",dpi=500,transparent=true)
    display(fig)
    close(fig)
    return nothing 
end

function plot_density_matrix_ivc_spin(metadata::String;ik::Int=1)
    # only works for (-2,-2)
    # plotting the spin for the two decoupled ivcs, namely {K↑,K'↓} versus {K↓,K'↑},
    # define spin operator as Sx τx
    hf = load(metadata,"hf");
    s0 = ComplexF64[1 0;0 1]
    sx = ComplexF64[0 1;1 0]
    sy = ComplexF64[0 -1im;1im 0]
    sz = ComplexF64[1 0;0 -1]
    Iq = Array{ComplexF64}(I,2hf.q,2hf.q)

    sxτxop = kron(sx,sx)
    syτxop = kron(sy,sx)
    sxτyop = kron(sx,sy)
    syτyop = kron(sy,sy)
    
    P = reshape(hf.P,2hf.q,4,2hf.q,4,:)[:,:,:,:,ik:ik]
    
    P1 = copy(P)
    P1[:,[2,3],:,[2,3],:] .= 0.0
    P2 = copy(P)
    P2[:,[1,4],:,[1,4],:] .= 0.0

    sxτx1 = reshape( sum(P1.*reshape(sxτxop,1,4,1,4,1),dims=(2,4,5)) / (4size(hf.P,3)), 2hf.q,2hf.q )
    sxτy1 = reshape( sum(P1.*reshape(sxτyop,1,4,1,4,1),dims=(2,4,5)) / (4size(hf.P,3)) , 2hf.q,2hf.q)
    syτx1 = reshape( sum(P1.*reshape(syτxop,1,4,1,4,1),dims=(2,4,5)) / (4size(hf.P,3)) , 2hf.q,2hf.q)
    syτy1 = reshape( sum(P1.*reshape(syτyop,1,4,1,4,1),dims=(2,4,5)) / (4size(hf.P,3)) , 2hf.q,2hf.q)

    sxτx2 = reshape( sum(P2.*reshape(sxτxop,1,4,1,4,1),dims=(2,4,5)) / (4size(hf.P,3)) , 2hf.q,2hf.q)
    sxτy2 = reshape( sum(P2.*reshape(sxτyop,1,4,1,4,1),dims=(2,4,5)) / (4size(hf.P,3)) , 2hf.q,2hf.q)
    syτx2 = reshape( sum(P2.*reshape(syτxop,1,4,1,4,1),dims=(2,4,5)) / (4size(hf.P,3)) , 2hf.q,2hf.q)
    syτy2 = reshape( sum(P2.*reshape(syτyop,1,4,1,4,1),dims=(2,4,5)) / (4size(hf.P,3)) , 2hf.q,2hf.q)

    return [sxτx1,sxτy1,syτx1,syτy1], [sxτx2,sxτy2,syτx2,syτy2]

end