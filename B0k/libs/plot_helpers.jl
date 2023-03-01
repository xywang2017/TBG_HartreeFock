function plot_energies(ϵk::Vector{Float64},chern::Vector{Float64};lines::Vector{Float64}=[])
    fig = figure(figsize=(3,3))
    pl = scatter(zeros(length(ϵk)),ϵk,s=6,c=chern,cmap="coolwarm",vmin=-1,vmax=1)
    colorbar(pl)
    for line in lines 
        axhline(line,ls=":",c="gray")
    end
    xlim([-0.1,0.1])
    ylabel("E (meV)")
    # ylim([-30,30])
    tight_layout()
    savefig("test.pdf")
    display(fig)
    close(fig)
    return nothing
end