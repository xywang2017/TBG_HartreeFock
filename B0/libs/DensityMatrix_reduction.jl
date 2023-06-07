include(joinpath(fpath,"B0/libs/HFChern_mod.jl"))

mutable struct DensityMat 
    P::Array{ComplexF64,3} # original density matrix
    δs::Array{Float64,2} # density matrix decomposition coefficients
    basis_mat::Vector{Matrix{ComplexF64}} # vector container of decomposition matrices 

    corr::Array{Float64,2} # correlation matrix defined as sum_{k} 
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

    norm_diff = norm(P1 .- A.P)
    println("Norm difference of reconstructed density matrix is: ",norm_diff)
    return norm_diff
end