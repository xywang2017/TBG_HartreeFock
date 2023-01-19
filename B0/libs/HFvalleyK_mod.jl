include("BM_mod.jl")

mutable struct HartreeFock
    bm::HBM # HBM 

    ν::Float64 # filling fraction -4 to 4 
    P::Array{ComplexF64,3}  # one particle density matrix
    ϵk::Matrix{Float64} # eigenvalues of HF renormalized band dispersions nfl x ns x nk

    H::Array{ComplexF64,3} # running Hamiltonian for any given k; nfl*ns x nfl*ns x nk
    precision::Float64  # iteration stopping point

    HartreeFock() = new()
end

@inline function V(q::ComplexF64,Lm::Float64) ::Float64
    res = 1e-6
    ϵr = 5.0
    # return  ( abs(q) < res ) ? 0 : 2π/(ϵr*abs(q))
    return ( abs(q) < res ) ? 0 : 2π/(ϵr*abs(q))*tanh(abs(q)*Lm/2)
end

function run_HartreeFock(hf::HartreeFock,bm::HBM,ν::Float64=0.0)
    hf.bm = bm 
    hf.ν = ν 
    hf.precision = 1e-8

    hf.P = zeros(ComplexF64,bm.nb,bm.nb,bm.latt.nk)
    hf.H = zeros(ComplexF64,bm.nb,bm.nb,bm.latt.nk)
    hf.ϵk = collect(reshape(bm.hbm,bm.nb,bm.latt.nk))

    init_P(hf.P,hf.ν)

    ## Hartree Fock contributions 
    Λkp = blk.Λkp

    kvec = real(bm.latt.kvec)*params.g1 + imag(bm.latt.kvec)*params.g2

    # Coulomb scale in meV units 
    ee = 1.6e-19
    ϵϵ = 8.8541878128e−12	
    aa = 2.46e-10
    area_moire = abs(bm.params.a1'*bm.params.a2)
    V0 = ee/(4π*ϵϵ*area_moire*aa) * 1e3
    
    ## H F iterations
    norm_convergence = 10.0 
    iter = 0
    H0 = zeros(ComplexF64,bm.nb,bm.nb,bm.latt.nk)
    add_BM(H0,bm.hbm,bm.nb,bm.latt.nk;α=1.0)
    iter_err = Float64[]
    iter_energy = Float64[]
    while norm_convergence > hf.precision
        hf.H .= H0
        add_Hartree(hf.H,hf.P,Λkp,bm.Gs,bm.nb,bm.latt.nk;β=1.0,V0=V0)
        add_Fock(hf.H,hf.P,Λkp,bm.Gs,kvec,bm.nb,bm.latt.nk;β=1.0,V0=V0)
        Etot = compute_HF_energy(hf.H .- H0,H0,hf.P)
        #Δ is a projector to make it closed shell
        if norm_convergence <5e-5
            Δ = 0.0 
        else 
            Δ = 0.01 #min(iter * 0.2,1)
        end
        norm_convergence = update_P(hf.P,hf.H,hf.ϵk,hf.ν;Δ=Δ)
        
        # tmpP = zeros(ComplexF64,bm.nfl*bm.ns,bm.nfl*bm.ns,bm.latt.nk)
        # if (iter%50==0)
        #     init_P(tmpP,hf.ν)
        #     hf.P .= 0.9*hf.P .+ 0.1*tmpP
        # end
        iter +=1
        println("Iter: ",iter)
        println("Running HF energy: ",Etot)
        println("Running norm convergence: ",norm_convergence)
        push!(iter_energy,compute_HF_energy(hf.H .- H0,H0,hf.P))
        push!(iter_err,norm_convergence)
        
    end

    return nothing
end

function init_P(P::Array{ComplexF64,3},ν::Float64; _RandInit::Bool=true,P0::Array{ComplexF64,3}=ones(ComplexF64,1,1,1))
    """
        Initialization of density matrix 
    """
    if _RandInit
        dim_to_keep = round(Int,size(P,1) * size(P,3) * (ν+1)/2)
        nk_to_populate = randperm(size(P,1) * size(P,3))[1:dim_to_keep]
        for ii in nk_to_populate
            ik = (ii-1)÷2 + 1
            n = (ii-1)%2 + 1
            P[n,n,ik] = 1 - 0.5
            P[3-n,3-n,ik] = 0 - 0.5
            # check_Hermitian(P[:,:,ik])
        end
    else
        P .= P0 
    end
    # println("Trace of initP is: ", sum([tr(P[:,:,ik]) for ik in 1:size(P,3)]) )
    return nothing
end

function add_BM(H::Array{ComplexF64,3},hbm::Array{Float64},nb::Int,nk::Int;α::Float64=1.0)
    """ 
        Add BM kinetic energy term suppressed by parameter α
    """
    
    tmp_hbm = reshape(hbm,nb,nk)
    tmp_H = reshape(H,nb,nb,nk)
    @inbounds @fastmath begin 
        for ifl in 1:nb
            tmp_H[ifl,ifl,:] .+= α * tmp_hbm[ifl,:]
        end
    end
    return nothing
end

function add_Hartree(H::Array{ComplexF64,3},P::Array{ComplexF64,3},Λkp::Array{ComplexF64,3},Gs::Vector{ComplexF64},
                       nb::Int,nk::Int;β::Float64=1.0,V0::Float64=1.0)
    """
        Hartree Contribution suppressed by parameter β
    """
    gMoire = abs(Gs[(length(Gs)+1)÷2+1])
    Lm = 4π/(sqrt(3)*gMoire)
    Λdiag = zeros(ComplexF64,nb,nb) 
    trPG = zeros(ComplexF64,length(Gs))
    tmp_Λkp = reshape(Λkp,nb,nk,nb,nk,length(Gs))
    @inbounds @fastmath begin 
        for ig in eachindex(Gs)
            for ik in 1:nk 
                Λdiag .= view(tmp_Λkp,:,ik,:,ik,ig)
                trPG[ig] += tr(view(P,:,:,ik)*conj(Λdiag)) *4 # spin-valley 4
            end
        end
        for ik in 1:size(H,3) 
            for ig in eachindex(Gs)
                H[:,:,ik] .+= ( β/nk * V0*V(Gs[ig],Lm) * trPG[ig]) * view(tmp_Λkp,:,ik,:,ik,ig)
            end
        end
    end
    return nothing
end

function add_Fock(H::Array{ComplexF64,3},P::Array{ComplexF64,3},Λkp::Array{ComplexF64,3},Gs::Vector{ComplexF64},kvec::Vector{ComplexF64},
                    nb::Int,nk::Int;β::Float64=1.0,V0::Float64=1.0)
    """
        Fock Contribution 
    """
    tmp_Λkp = reshape(Λkp,nb,nk,nb,nk,length(Gs))
    tmp_Fock = zeros(ComplexF64,nb,nb)
    gMoire = abs(Gs[(length(Gs)+1)÷2+1])
    Lm = 4π/(sqrt(3)*gMoire)
    @inbounds @fastmath begin 
        for ik in 1:nk
            tmp_Fock .= 0.0 + 0.0im
            for ig in eachindex(Gs), ip in 1:nk
                tmp_Fock .+= ( β*V0*V(kvec[ip]-kvec[ik]+Gs[ig],Lm) /nk) * 
                            ( view(tmp_Λkp,:,ik,:,ip,ig)*transpose(view(P,:,:,ip))*view(tmp_Λkp,:,ip,:,ik,length(Gs)-ig+1) )
            end
            H[:,:,ik] .-= tmp_Fock
        end
    end
    return nothing
end

function add_Projector(H::Array{ComplexF64,3},P::Array{ComplexF64,3};Δ::Float64=10.0)
    """ 
        Add a projector to make it closed shell; Δ = 10meV
    """
    H .-= Δ .* P
    return nothing 
end

function update_P(P::Array{ComplexF64,3},H::Array{ComplexF64,3},ϵk::Matrix{Float64},ν::Float64;Δ::Float64=0.0)
    """
        Diagonalize Hamiltonian for every k; use ν to keep the lowest N particle states;
        update P 
    """
    @inbounds @fastmath begin
        νnorm = round(Int,(ν+1)/2 * size(H,1)*size(H,3))  # total number of occupied states 
        vecs = zeros(ComplexF64,size(H,1),size(H,2),size(H,3))
        for ik in 1:size(H,3)
            # check_Hermitian(H[:,:,ik])
            ϵk[:,ik],vecs[:,:,ik] = eigen(Hermitian(view(H,:,:,ik)-Δ*(conj.(view(P,:,:,ik))+0.5*I)) )
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

        # norm_convergence = norm(0.9*P_new .- 0.9*P) ./ norm(0.9*P_new+0.1*P)
        # P .= 0.9*P_new .+ 0.1*P
        norm_convergence = norm(P_new .- P) ./ norm(P_new)
        P .= P_new 
    end
    # println("Trace of P is: ", sum([tr(P_new[:,:,ik]) for ik in 1:size(P,3)]) )
    return norm_convergence
end

function compute_HF_energy(H_HF::Array{ComplexF64,3},H0::Array{ComplexF64,3},P::Array{ComplexF64,3})
    Etot = 0.0 
    for ik in 1:size(H0,3) 
        Etot += tr(view(H_HF,:,:,ik)*transpose(view(P,:,:,ik)))/2 + tr(view(H0,:,:,ik)*transpose(view(P,:,:,ik)))
    end
    return real(Etot)/(size(H0,3))
end

function check_Hermitian(H::Matrix{ComplexF64})
    err = norm(H - H')
    if err > 1e-6
        println("Error with Hermitian Hamiltonian")
    end
    return nothing
end