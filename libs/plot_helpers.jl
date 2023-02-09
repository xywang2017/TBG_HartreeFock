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

    pl=scatter(ones(length(ϵsorted))*0.25,ϵsorted,c=chern,cmap="Spectral_r",s=6,vmin=-1,vmax=1)
    colorbar(pl)
    ϵF = minimum(ϵsorted)
    δν = 1/length(ϵsorted) * 8
    Δ = 0.0 
    for i in 2:length(ϵsorted)
        ν = (i / length(ϵsorted) -0.5) * 8 
        if ν>=νF+δν/2 
            ϵF = (ν==νF) ? (ϵsorted[i] + ϵsorted[i+1])/2 : (ϵsorted[i] + ϵsorted[i-1])/2
            Δ = (ν==νF) ? (ϵsorted[i+1] - ϵsorted[i]) : (ϵsorted[i] - ϵsorted[i-1])
            println("Gap size: ", Δ)
            break 
        end
    end
    axhline(ϵF,ls=":",c="gray")
    ylabel("E (meV)")
    xlabel(L"ϕ/ϕ_0")
    # legend()
    # ylim([-0.4,0.8])
    tight_layout()
    # savefig(savename,dpi=500)
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