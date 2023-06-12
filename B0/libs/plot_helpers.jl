function plot_contour_maps(kvec::Matrix{ComplexF64},ϵ::Matrix{Float64};
        points::Vector{ComplexF64}=[],contourlines::Vector{Float64}=[],limits::Vector{Float64}=Float64[])
    kx,ky = real(kvec), imag(kvec)
    fig = figure(figsize=(4,3))
    if !isempty(limits)
        pl=contourf(kx,ky,ϵ,cmap="bwr",levels=20,vmin=limits[1],vmax=limits[2])
    else
        # bound = maximum(abs.(ϵ))
        pl=contourf(kx,ky,ϵ,cmap="Spectral_r",levels=20)
    end
    if !isempty(contourlines)
        contour(kx,ky,ϵ,levels=contourlines)
    end
    plot(real(points),imag(points),"k+")
    colorbar(pl)
    xlabel(L"k_1")
    ylabel(L"k_2")
    axis("equal")
    tight_layout()
    savefig("test.pdf",transparent=true)
    display(fig)
    close(fig)
    return nothing
end

function plot_density_maps(kvec::Matrix{ComplexF64},ϵ::Matrix{Float64};
    points::Vector{ComplexF64}=[],contourlines::Vector{Float64}=[],limits::Vector{Float64}=Float64[])
    kx,ky = real(kvec), imag(kvec)
    fig = figure(figsize=(4,3))
    if !isempty(limits)
        pl=pcolormesh(kx,ky,ϵ,cmap="bwr",vmin=limits[1],vmax=limits[2])
    else
        bound = maximum(abs.(ϵ))
        pl=pcolormesh(kx,ky,ϵ,cmap="bwr",vmin=-bound,vmax=bound)
    end
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


function plot_density_maps_collective(kvec::Matrix{ComplexF64},Δs::Array{Float64,5};
    points::Vector{ComplexF64}=[],contourlines::Vector{Float64}=[],limits::Vector{Float64}=Float64[])
    kx,ky = real(kvec), imag(kvec)
    γstr = ["γx","γy","γz"]
    τstr = ["τx","τy","τz"]
    sstr = ["sx","sy","sz"]
    titlestrs = [γstr,τstr,sstr]
    fig,ax = subplots(3,3,sharex=true,sharey=true,figsize=(12,9))
    for r in 1:3, c in 1:3
        str = titlestrs[c][r]
        if c == 1
            ϵ = Δs[r,4,4,:,:]
        elseif c == 2
            ϵ = Δs[4,r,4,:,:]
        elseif c ==3 
            ϵ = Δs[4,4,r,:,:]
        end
        bound = maximum(abs.(ϵ))
        pl=ax[r,c].pcolormesh(kx,ky,ϵ,cmap="bwr",vmin=-bound,vmax=bound)
        # pl=ax[c,r].pcolormesh(kx,ky,ϵ,cmap="bwr")
        ax[r,c].set_title(str)
        colorbar(pl,ax=ax[r,c])
    end
    for r in 1:3, c in 1:3
        if r == 3 
            ax[r,c].set_xlabel(L"k_x/|g_1|")
        end
        if c == 1
            ax[r,c].set_ylabel(L"k_y/|g_1|")
        end
    end
    # ax[1,1].set_aspect("equal", adjustable="datalim")
    savefig("test.pdf")
    display(fig)
    close(fig)
    return nothing
end


function plot_density_maps_collectivev0(kvec::Matrix{ComplexF64},Δs::Array{Float64,5};
    points::Vector{ComplexF64}=[],contourlines::Vector{Float64}=[],limits::Vector{Float64}=Float64[])
    kx,ky = real(kvec), imag(kvec)
    γstr = ["γ0","γx","γy","γz"]
    τstr = ["τ0","τx","τy","τz"]
    sstr = ["s0","sx","sy","sz"]
    fig,ax = subplots(4,4,sharex=true,sharey=true,figsize=(16,16))
    order_params = view(Δs,:,:,1,:,:)
    bound = maximum(abs.(order_params))
    for r in 1:4, c in 1:4
        ϵ = order_params[r,c,:,:]
        bound = maximum(abs.(ϵ))
        str = γstr[r]*τstr[c]
        # str = τstr[r]*sstr[c]
        pl=ax[r,c].pcolormesh(kx,ky,ϵ,cmap="bwr",vmin=-bound,vmax=bound)
        # pl=ax[c,r].pcolormesh(kx,ky,ϵ,cmap="bwr")
        ax[r,c].set_title(str)
        colorbar(pl,ax=ax[r,c])
    end
    savefig("test.pdf")
    display(fig)
    close(fig)
    return nothing
end


function plot_energy_cutsv0(kvec::Vector{Float64},ϵ::Array{Float64,2};lines::Vector{Float64}=[])
    fig = figure(figsize=(4,3))
    colors = ["r","b"]
    valleys = ["K","K\'"]
    for i in 1:size(ϵ,1)
        if mod(i-1,2)==0
            plot(kvec,ϵ[i,:],"-",c=colors[(i-1)÷2+1],label=valleys[(i-1)÷2+1],lw=1)
        else
            plot(kvec,ϵ[i,:],"-",c=colors[(i-1)÷2+1],lw=1)
        end
    end
    for line in lines 
        axhline(line,ls=":",c="gray")
    end
    # xlim([minimum(kvec)-0.1,2.5*maximum(kvec)])
    xlabel("k")
    ylabel("E (meV)")
    # ylim([-30,30])
    legend()
    tight_layout()
    savefig("test.pdf",transparent=true)
    display(fig)
    close(fig)
    return nothing
end



function plot_energy_cuts(kvec::Vector{Float64},ϵ::Array{Float64,2};lines::Vector{Float64}=[])
    fig = figure(figsize=(3,3))
    for i in 1:size(ϵ,1)
        plot(kvec,ϵ[i,:],"-",lw=1)
    end
    for line in lines 
        axhline(line,ls=":",c="gray")
    end
    # xlim([minimum(kvec)-0.1,2.5*maximum(kvec)])
    xlabel("k")
    ylabel("E (meV)")
    # ylim([-30,30])
    tight_layout()
    savefig("test.pdf",transparent=true)
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

