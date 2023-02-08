function plot_contour_maps(kvec::Matrix{ComplexF64},ϵ::Matrix{Float64};
        points::Vector{ComplexF64}=[])
    kx,ky = real(kvec), imag(kvec)
    fig = figure(figsize=(4,3))
    pl=contourf(kx,ky,ϵ,cmap="Spectral_r",levels=20)
    plot(real(points),imag(points),"k+")
    colorbar(pl)
    xlabel(L"k_x")
    ylabel(L"k_y")
    axis("equal")
    tight_layout()
    display(fig)
    close(fig)
    return nothing
end

function plot_energy_cuts(kvec::Vector{Float64},ϵ::Array{Float64,2})
    fig = figure(figsize=(4,3))
    for i in 1:size(ϵ,1)
        plot(kvec,ϵ[i,:],"o-",ms=2,markeredgecolor="none",label="band $(i)")
    end
    xlabel("k")
    ylabel("E (meV)")
    legend()
    tight_layout()
    display(fig)
    close(fig)
    return nothing
end