include("BM_mod.jl")
using Tullio

mutable struct HartreeFock
    bm::HBM # HBM 

    ν::Float64 # filling fraction -4 to 4 
    P::Array{ComplexF64,3}  # one particle density matrix
    ϵk::Matrix{Float64} # eigenvalues of HF renormalized band dispersions nfl x ns x nk

    H::Array{ComplexF64,3} # running Hamiltonian for any given k; nfl*ns x nfl*ns x nk
    precision::Float64  # iteration stopping point

    HartreeFock() = new()
end

@inline function V(q::ComplexF64) ::Float64
    res = 1e-6
    ϵr = 9.5
    return  ( abs(q) < res ) ? 0 : 2π/(ϵr*abs(q))
end

function run_HartreeFock(hf::HartreeFock,bm::HBM,ν::Float64=0.0)
    hf.bm = bm 
    hf.ν = ν 
    hf.precision = 1e-5

    hf.P = zeros(ComplexF64,bm.nfl*bm.ns,bm.nfl*bm.ns,bm.latt.nk)
    hf.H = zeros(ComplexF64,bm.nfl*bm.ns,bm.nfl*bm.ns,bm.latt.nk)
    hf.ϵk = zeros(Float64,bm.nfl*bm.ns,bm.latt.nk)
    tmp_ϵk = reshape(hf.ϵk,bm.ns,bm.nfl*bm.latt.nk)
    tmp_ϵk[1,:] = bm.hbm 
    tmp_ϵk[2,:] = bm.hbm

    init_P(hf.P,hf.ν)

    ## Hartree Fock contributions 
    Λkp = zeros(ComplexF64,bm.ns*bm.nfl*bm.latt.nk,bm.ns*bm.nfl*bm.latt.nk,length(bm.Gs))
    tmp_Λkp = reshape(Λkp,bm.ns,bm.nfl*bm.latt.nk,bm.ns,bm.nfl*bm.latt.nk,length(bm.Gs))
    tmp_Λkp[1,:,1,:,:] = blk.Λkp 
    tmp_Λkp[2,:,2,:,:] = blk.Λkp

    kvec = real(bm.latt.kvec)*params.g1 + imag(bm.latt.kvec)*params.g2

    ## H F iterations
    norm_convergence = 10 
    iter = 0
    while norm_convergence > hf.precision
        hf.H .= 0.0 + 0.0im
        add_BM(hf.H,bm.hbm,bm.ns,bm.nη,bm.nb,bm.latt.nk;α=0.0)
        add_Hartree(hf.H,hf.P,Λkp,bm.Gs,bm.ns,bm.nη,bm.nb,bm.latt.nk;β=1.0)
        add_Fock(hf.H,hf.P,Λkp,bm.Gs,kvec,bm.ns,bm.nη,bm.nb,bm.latt.nk;β=1.0)

        norm_convergence = update_P(hf.P,hf.H,hf.ϵk,hf.ν)
        iter +=1
        println("Iter: ",iter," Running norm convergence: ",norm_convergence)
    end

    return nothing
end

function init_P(P::Array{ComplexF64,3},ν::Float64; _RandInit::Bool=true,P0::Array{ComplexF64,3}=ones(ComplexF64,1,1,1))
    """
        Initialization of density matrix 
    """
    if _RandInit
        mat = zeros(ComplexF64,size(P,1),size(P,2))
        dim_to_keep = round(Int,size(P,1) * size(P,3) * (ν+4)/8)
        for ik in 1:size(P,3)
            # mat .= rand(ComplexF64,size(P,1),size(P,2))
            mat .= kron(rand(ComplexF64,size(P,1)÷2,size(P,2)÷2),ComplexF64[1 0;0 1])  # assuming spin degeneracy
            vec = eigvecs(Hermitian(mat))[:,1:(dim_to_keep÷size(P,3))]
            P[:,:,ik] = conj(vec) * transpose(vec) - 0.5*I
            # check_Hermitian(P[:,:,ik])
        end
    else
        P .= P0 
    end
    # println("Trace of initP is: ", sum([tr(P[:,:,ik]) for ik in 1:size(P,3)]) )
    return nothing
end

function add_BM(H::Array{ComplexF64,3},hbm::Array{Float64},
                    ns::Int,nη::Int,nb::Int,nk::Int;α::Float64=1.0)
    """ 
        Add BM kinetic energy term suppressed by parameter α
    """
    tmp_hbm = reshape(hbm,nη*nb,nk)
    tmp_H = reshape(H,ns,nη*nb,ns,nη*nb,nk)
    for ifl in 1:nη*nb
        tmp_H[1,ifl,1,ifl,:] .+= α * tmp_hbm[ifl,:]
        tmp_H[2,ifl,2,ifl,:] .+= α * tmp_hbm[ifl,:]
    end
    return nothing
end

function add_Hartree(H::Array{ComplexF64,3},P::Array{ComplexF64,3},Λkp::Array{ComplexF64,3},Gs::Vector{ComplexF64},
                        ns::Int,nη::Int,nb::Int,nk::Int;β::Float64=1.0)
    """
        Hartree Contribution suppressed by parameter β
    """
    trPG = zeros(ComplexF64,length(Gs))
    tmp_Λkp = reshape(Λkp,ns*nη*nb,nk,ns*nη*nb,nk,length(Gs))
    for ig in eachindex(Gs)
        for ik in 1:nk 
            trPG[ig] += tr(view(P,:,:,ik)*conj(view(tmp_Λkp,:,ik,:,ik,ig)))
        end
    end
    for ik in 1:size(H,3) 
        for ig in eachindex(Gs)
            H[:,:,ik] .+= ( β/nk * V(Gs[ig]) * trPG[ig]) * view(tmp_Λkp,:,ik,:,ik,ig)
        end
    end
    return nothing
end

function add_Fock(H::Array{ComplexF64,3},P::Array{ComplexF64,3},Λkp::Array{ComplexF64,3},Gs::Vector{ComplexF64},kvec::Vector{ComplexF64},
                        ns::Int,nη::Int,nb::Int,nk::Int;β::Float64=1.0)
    """
        Fock Contribution 
    """
    tmp_Λkp = reshape(Λkp,ns*nη*nb,nk,ns*nη*nb,nk,length(Gs))
    for ik in 1:nk
        for ig in eachindex(Gs)
            for ip in 1:nk
                H[:,:,ik] .-= ( β*V(kvec[ip]-kvec[ik]+Gs[ig]) /nk) * 
                            ( view(tmp_Λkp,:,ik,:,ip,ig)*transpose(view(P,:,:,ip))*view(tmp_Λkp,:,ip,:,ik,length(Gs)-ig+1) )
            end
        end
    end
    # ntot = ns*nη*nb
    # Vkp_mat = reshape((β/nk).* V.(reshape(kvec,1,:,1).-reshape(kvec,:,1,1) .+ reshape(Gs,1,1,:) ),
    #                     1,1,nk,1,1,nk,length(Gs))
    # λl = reshape(Λkp,ntot,1,nk,1,ntot,nk,length(Gs))
    # Pvec = reshape(P,1,1,1,ntot,ntot,nk,1)
    # λr = reshape(Λkp,1,ntot,nk,ntot,1,nk,length(Gs))
    # tmp_FockG = zeros(ComplexF64,ntot,ntot,nk,ntot,ntot,nk)
    # for ig in eachindex(Gs)
    #     tmp_FockG .= view(Vkp_mat,:,:,:,:,:,:,ig) .*view(λl,:,:,:,:,:,:,ig) .* 
    #                  view(Pvec,:,:,:,:,:,:,1) .* conj.(view(λr,:,:,:,:,:,:,ig))
    #     H .-= sum(tmp_FockG,dims=(4,5,6))
    # end

    # @tullio H[α,β,ik] += -Vkp_mat[ik,ip,iG] * λkp[α,ik,γ1,ip,iG] * P[γ2,γ1,ip] * conj(λkp[β,ik,γ2,ip,iG])
    
    return nothing
end

function update_P(P::Array{ComplexF64,3},H::Array{ComplexF64,3},ϵk::Matrix{Float64},ν::Float64)
    """
        Diagonalize Hamiltonian for every k; use ν to keep the lowest N particle states;
        update P 
    """
    νnorm = round(Int,(ν+4)/8 * size(H,1)*size(H,3))  # total number of occupied states 

    vecs = zeros(ComplexF64,size(H,1),size(H,2),size(H,3))
    for ik in 1:size(H,3)
        # check_Hermitian(H[:,:,ik])
        ϵk[:,ik],vecs[:,:,ik] = eigen(Hermitian(view(H,:,:,ik)))
    end
    iϵ_sorted = sortperm(ϵk[:])
    iϵ_occupied = iϵ_sorted[1:νnorm]
    iband_occupied = (iϵ_occupied .-1) .% size(ϵk,1) .+1
    ik_occupied = (iϵ_occupied .-1) .÷ size(ϵk,1) .+1

    P_new = zeros(ComplexF64,size(P,1),size(P,2),size(P,3))
    for ik in 1:size(P,3)
        occupied_vecs = vecs[:,iband_occupied[ik_occupied.==ik],ik]
        P_new[:,:,ik] = conj(occupied_vecs)*transpose(occupied_vecs) - 0.5*I
    end

    norm_convergence = [norm(P_new[:,:,ik] .- P[:,:,ik]) ./ norm(P_new[:,:,ik]) for ik in 1:size(P,3)]
    P .= P_new
    # println("Trace of P is: ", sum([tr(P_new[:,:,ik]) for ik in 1:size(P,3)]) )
    return sum(norm_convergence)/length(norm_convergence)
end

function check_Hermitian(H::Matrix{ComplexF64})
    err = norm(H - H')
    if err > 1e-6
        println("Error with Hermitian Hamiltonian")
    end
    return nothing
end