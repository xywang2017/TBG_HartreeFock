using PyPlot 

## plot Hartree Fock spectra
function plot_spectra(ϵk::Matrix{Float64},σzτz::Matrix{Float64},νF::Float64,params::Params;savename::String="tmp.pdf")
    ee = 1.6e-19
    ϵϵ = 8.8541878128e−12	
    aa = 2.46e-10
    ϵr = 5.0
    Vcoulomb = ee/(4π*ϵϵ*ϵr* abs(params.a1)*aa) * 1e3
    
    fig = figure(figsize=(3,3))
    idx = sortperm(ϵk[:])
    ϵsorted = ϵk[idx] #./Vcoulomb
    chern = σzτz[idx]

    # ϵsorted = ϵsorted[chern .>0] #./Vcoulomb
    # chern = chern[chern .>0]

    pl=scatter(ones(length(ϵsorted))*0.25,ϵsorted,c=chern,cmap="coolwarm",s=6,vmin=-1,vmax=1)
    colorbar(pl)
    ν = eachindex(ϵsorted) ./ length(ϵsorted)
    i = 1
    while (νF+4)/8 > ν[i]
        i += 1
    end
    ϵF = (ϵsorted[i+1] + ϵsorted[i])/2 
    Δ = (ϵsorted[i+1] - ϵsorted[i]) 
    
    println("Gap size: ", Δ)
    axhline(ϵF,ls=":",c="gray")
    ylabel("E (meV)")
    xlabel(L"ϕ/ϕ_0")
    # legend()
    # ylim([-0.4,0.8])
    tight_layout()
    savefig(savename,transparent=true)
    display(fig)
    close(fig)
    println("Sublattice polarization operator is: ",sum(chern[ϵsorted.<ϵF])/(size(ϵk,2)*size(ϵk,1))*8)
    # println("Total energy: ",(0.25*sum(ϵsorted[ϵsorted.<ϵF])-0.25*sum(ϵsorted[ϵsorted.>=ϵF]))/size(hf.ϵk,2)/size(hf.ϵk,1)*8)
    # return sum(chern[ϵsorted.<ϵF])/(size(ϵk,2)*size(ϵk,1))*8
    return Δ /Vcoulomb
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
function plot_density_matrix_bm(fname::String)
    # plot at a given k point, 2qx4 x 2qx4 matrix 
    hf = load(fname,"hf");
    ik = 1
    fig = figure(figsize=(6,6))
    P0 = view(hf.P,:,:,ik) + 0.5I
    pl = imshow(abs.(P0),vmin=0,vmax=1)
    colorbar(pl)
    tight_layout()
    display(fig)
    close(fig)
    return nothing
end


### density matrix analysis 
function plot_density_matrix_bm_valley_spin(fname::String)
    # plot at a given k point, 2qx4 x 2qx4 matrix 
    hf = load(fname,"hf");
    ik = 1
    fig = figure(figsize=(6,6))
    P0 = reshape(view(hf.P,:,:,ik)+0.5I,2hf.q,4,2hf.q,4);
    fig,ax = subplots(2,2,figsize=(8,8))
    states = ["K↑","K'↑","K↓","K'↓"]
    for r in 1:2, c in 1:2 
        pl=ax[r,c].imshow(abs.(P0[:,r+2(c-1),:,r+2(c-1)]),vmin=0,vmax=1)
        colorbar(pl,ax=ax[r,c],fraction=0.046, pad=0.04)
        ax[r,c].set_title(states[r+2(c-1)])
    end
    tight_layout()
    display(fig)
    close(fig)
    return nothing
end

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
    pl=imshow(abs.(Pstrong),vmin=0,vmax=1)
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
    fig,ax = subplots(2,2,figsize=(8,8))
    states = ["K↑","K'↑","K↓","K'↓"]
    for r in 1:2, c in 1:2 
        pl=ax[r,c].imshow(abs.(Pstrong[:,:,r+2(c-1)]),vmin=0,vmax=1)
        colorbar(pl,ax=ax[r,c],fraction=0.046, pad=0.04)
        ax[r,c].set_title(states[r+2(c-1)])
    end
    tight_layout()
    display(fig)
    close(fig)
    return nothing 
end