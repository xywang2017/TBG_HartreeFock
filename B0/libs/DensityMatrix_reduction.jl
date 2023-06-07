include(joinpath(fpath,"B0/libs/HFChern_mod.jl"))

mutable struct DensityMat 
    P::Array{ComplexF64,3} # original density matrix
    δs::Array{Float64,2} # density matrix decomposition coefficients
    basis_mat::Vector{Matrix{ComplexF64}} # vector container of decomposition matrices 

    corr::Array{Float64,2} # correlation matrix defined as sum_{k} 
    corr_vals::Vector{Float64} # correlation matrix eigenvalues
    φs::Array{Float64,2} # extracted basis vectors from analysis of correlation matrix 
    U::Array{Float64,2} # conversion coefficients defined as δs_a = sum_{b} U_{a,b}φs_{b}
    Oφs::Array{ComplexF64,3} # 8x8 matrices associated with each ϕs

    DensityMat() = new()
end

function constructDensityMat(hf::HartreeFock)
    A = DensityMat() 
    A.P = deepcopy(hf.P)
    
    s0 = ComplexF64[1 0;0 1]
    sx = ComplexF64[0 1;1 0]
    sy = ComplexF64[0 -1im;1im 0]
    sz = ComplexF64[1 0;0 -1]
    paulis = [s0,sx,sy,sz]
    A.basis_mat = []
    for i in 1:4, j in 1:4, k in 1:4
        push!(A.basis_mat, kron(paulis[i],kron(paulis[j],paulis[k])) )
    end

    A.δs = zeros(Float64,64,size(A.P,3))
    for ik in 1:size(A.P,3),j in 1:64
        A.δs[j,ik] = real(tr(view(A.P,:,:,ik)*A.basis_mat[j])/8 )
    end

    # --------------------- reconstruction --------------------- # 
    A.corr = ( A.δs * A.δs' )/size(A.δs,2)
    F = eigen(A.corr)
    A.U = F.vectors 
    A.corr_vals = F.values
    A.φs = F.vectors'* A.δs
    A.Oφs = zeros(ComplexF64,8,8,64)
    for α in 1:64, j in 1:64
        A.Oφs[:,:,α] += F.vectors[j,α]*A.basis_mat[j]
    end

    # checkReconstructionValidity(A,collect(50:64))
    return A
end

function checkReconstructionValidity(A::DensityMat,reconstr_range::Vector{Int})
    P1 = zeros(ComplexF64,8,8,size(A.P,3))
    for ik in 1:size(A.P,3)
        for α in reconstr_range
            P1[:,:,ik] .+= A.φs[α,ik] * view(A.Oφs,:,:,α)
        end
    end

    norm_diff = norm(P1 .- A.P)/size(A.P,3)
    println("Norm difference of reconstructed density matrix is: ",norm_diff)
    return norm_diff
end

function plot_formfactor_info(dm::DensityMat,idx::Int)
    fig,ax = subplots(3,1,figsize=(4,9))
    pl = ax[1].pcolormesh(real(kvec),imag(kvec),reshape(dm.φs[idx,:],lk,lk),cmap="bwr")
    colorbar(pl,ax=ax[1])
    ax[1].set_xlabel(L"k_1")
    ax[1].set_ylabel(L"k_2")
    ax[1].set_title("φ(k)")
    ax[1].axis("equal")
    bound = maximum(abs.(dm.Oφs[:,:,idx]))
    pl = ax[2].imshow(real.(dm.Oφs[:,:,idx]),origin="lower",extent=(0.5,8.5,0.5,8.5),vmin=-bound,vmax=bound,cmap="bwr")
    ax[2].axhline(4.5,c="k",ls=":")
    ax[2].axvline(4.5,c="k",ls=":")
    ax[2].set_title("Re[Oφ]")
    colorbar(pl,ax=ax[2])
    pl = ax[3].imshow(imag.(dm.Oφs[:,:,idx]),origin="lower",extent=(0.5,8.5,0.5,8.5),vmin=-bound,vmax=bound,cmap="bwr")
    ax[3].set_title("Im[Oφ]")
    ax[3].axhline(4.5,c="k",ls=":")
    ax[3].axvline(4.5,c="k",ls=":")
    colorbar(pl,ax=ax[3])
    tight_layout()
    savefig("test.pdf",transparent=true)
    display(fig)
    close(fig)
    return nothing
end



function plot_corr_values(dm::DensityMat)
    fig = figure(figsize=(4,3))
    plot(eachindex(dm.corr_vals),dm.corr_vals,"g.")
    yscale("log")
    tight_layout()
    display(fig)
    close(fig)
    return nothing
end