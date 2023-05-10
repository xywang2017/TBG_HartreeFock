function plot_contour_maps(kvec::Matrix{ComplexF64},ϵ::Matrix{Float64};
        points::Vector{ComplexF64}=[],contourlines::Vector{Float64}=[])
    kx,ky = real(kvec), imag(kvec)
    fig = figure(figsize=(4,3))
    pl=contourf(kx,ky,ϵ,cmap="Spectral_r",levels=20)
    if !isempty(contourlines)
        contour(kx,ky,ϵ,levels=contourlines)
    end
    plot(real(points),imag(points),"k+")
    colorbar(pl)
    xlabel(L"k_x/|g_1|")
    ylabel(L"k_y/|g_1|")
    axis("equal")
    tight_layout()
    savefig("test.pdf")
    display(fig)
    close(fig)
    return nothing
end

function plot_energy_cuts(kvec::Vector{Float64},ϵ::Array{Float64,2};lines::Vector{Float64}=[])
    fig = figure(figsize=(4,4))
    for i in 1:size(ϵ,1)
        plot(kvec,ϵ[i,:],"-",ms=2,markeredgecolor="none",label="band $(i)")
    end
    for line in lines 
        axhline(line,ls=":",c="gray")
    end
    # xlim([minimum(kvec)-0.1,2.5*maximum(kvec)])
    xlabel("k")
    ylabel("E (meV)")
    # ylim([-30,30])
    # legend()
    tight_layout()
    savefig("test.pdf")
    display(fig)
    close(fig)
    return nothing
end

function plot_energy_cuts_with_order_parameters(kvec::Vector{Float64},ϵ::Array{Float64,2},σz::Array{Float64,2};lines::Vector{Float64}=[])
    fig = figure(figsize=(4,4))
    # for i in 1:size(ϵ,1)
    #     plot(kvec,ϵ[i,:],"-",c="gray",lw=1)
    # end
    pl = 0.0
    # for i in 1:size(ϵ,1)
    #     plot(kvec,ϵ[i,:],"-",c="gray")
    # end
    for i in 1:size(ϵ,1)
        pl = scatter(kvec,ϵ[i,:],s=6,c=σz[i,:],cmap="coolwarm",vmin=-1,vmax=1)
    end
    colorbar(pl)
    for line in lines 
        axhline(line,ls=":",c="gray")
    end
    xlim([minimum(kvec)-0.1,2.5*maximum(kvec)])
    xlabel("k")
    ylabel("E (meV)")
    # ylim([-30,30])
    # legend()
    tight_layout()
    savefig("test.pdf")
    display(fig)
    close(fig)
    return nothing
end

