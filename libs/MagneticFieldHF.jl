include("bmLL.jl")
using LinearAlgebra
using Random

mutable struct HartreeFock
    metadata::Vector{String}  # data files for BM related
    # 
    latt::Lattice 
    params::Params

    # flux related 
    p::Int 
    q::Int
    nq::Int # number of points in a moire zone 
    ng::Int # -3 to 3 

    # flavors 
    nb::Int 
    nη::Int 
    ns::Int 
    nfl::Int

    # gvec and overlap matrix for a given gvec
    gvec::Matrix{ComplexF64}
    Λ::Array{ComplexF64,3}

    # Coulomb unit 
    V0::Float64

    ν::Float64 # filling fraction -4 to 4 
    P::Array{ComplexF64,3}  # one particle density matrix
    ϵk::Matrix{Float64} # eigenvalues of HF renormalized band dispersions nfl x ns x nk
    σzτz::Matrix{Float64} # eigenvalues of σzτz for all HF states;

    Σz0::Array{ComplexF64,3} # projected Σz operator in the BM basis
    H0::Array{ComplexF64,3}  # non-interacting part of the Hamiltonian (BM)
    H::Array{ComplexF64,3} # running Hamiltonian for any given k; nfl*q x nfl*q x nk
    precision::Float64  # iteration stopping point

    HartreeFock() = new()
end

@inline function V(q::ComplexF64,Lm::Float64) ::Float64
    res = 1e-6
    ϵr = 10.0
    # return  ( abs(q) < res ) ? 0 : 2π/(ϵr*abs(q))
    return ( abs(q) < res ) ? 0 : 2π/(ϵr*abs(q))*tanh(abs(q)*Lm/2)
end

function CoulombUnit(params::Params)
     # Coulomb scale in meV units 
     ee = 1.6e-19
     ϵϵ = 8.8541878128e-12	
     aa = 2.46e-10
     area_moire = abs(imag(params.a1'*params.a2))
     V0 = ee/(4π*ϵϵ*area_moire*aa) * 1e3
     return V0
end

function run_HartreeFock(hf::HartreeFock,params::Params;precision::Float64=1e-5,
        ν::Float64=0.0,ϕ::Rational{Int}=1//10,prefix::String="",_Init::String="CNP",
        P0::Array{ComplexF64,3}=ones(ComplexF64,1,1,1),H0::Array{ComplexF64,3}=ones(ComplexF64,1,1,1))
    p, q = numerator(ϕ), denominator(ϕ)
    hf.p = p
    hf.q = q
    hf.ν = ν 
    hf.precision = precision
    hf.nb, hf.nη, hf.ns, hf.nfl = 2, 2, 2, 8
    hf.ng = 3
    hf.nq = (q>4) ? 2 : 2
    hf.metadata = [prefix*"_$(p)_$(q)_K_metadata.jld2",
                   prefix*"_$(p)_$(q)_Kprime_metadata.jld2"]
    
    hf.params = params
    hf.latt = Lattice() 
    constructLattice(hf.latt,hf.params;lk = hf.nq*hf.q)  #[0,1)x[0,1), so far works for p/q < 1
    hf.V0 = CoulombUnit(hf.params)

    # cannot save the overlap matrix for all of gvecs, so instead only save one gvec at a time
    hf.gvec = reshape(collect(-hf.ng:hf.ng),:,1)*hf.params.g1 .+ reshape(collect(-hf.ng*hf.q:hf.ng*hf.q),1,:)*hf.params.g2 ./hf.q
    hf.Λ = zeros(ComplexF64,2hf.q*hf.q*hf.nq^2,2hf.q*hf.q*hf.nq^2,hf.nη)

    # BM band structure and projected sublattice info
    hf.H0 = zeros(ComplexF64,hf.nfl*hf.q,hf.nfl*hf.q,hf.nq^2*hf.q)
    hf.Σz0 = zeros(ComplexF64,hf.nfl*hf.q,hf.nfl*hf.q,hf.nq^2*hf.q)
    BM_info(hf)

    # ------------------------------------------------- Begin Hartree Fock Procedure -------------------------------- #
    
    hf.P = zeros(ComplexF64,hf.nfl*hf.q,hf.nfl*hf.q,hf.nq^2*hf.q)
    hf.H = zeros(ComplexF64,hf.nfl*hf.q,hf.nfl*hf.q,hf.nq^2*hf.q)
    hf.ϵk = zeros(Float64,hf.nfl*hf.q,hf.nq^2*hf.q)
    hf.σzτz = zeros(Float64,hf.nfl*hf.q,hf.nq^2*hf.q)
    
    # --------- Initialization ---------- #
    
    init_P(hf,_Init=_Init,P0=P0,H0=H0)
    
    # --------- Hartree Fock Iterations ---------- #
    norm_convergence = 10.0 
    iter = 1
    iter_err = Float64[]
    iter_energy = Float64[]
    while norm_convergence > hf.precision
        println("Iter: ",iter)
        α = 1.0
        hf.H .= hf.H0 * α
        add_Hartree(hf;β=1.0,V0=hf.V0)
        # add_Fock(hf;β=1.0,V0=hf.V0)
        add_Fock_vectorize(hf;β=1.0,V0=hf.V0)
        Etot = compute_HF_energy(hf.H .- α*hf.H0,α*hf.H0,hf.P,hf.ν)
        #Δ is a projector to make it closed shell
        if norm_convergence <1e-4
            Δ = 0.0
        else 
            Δ = 0.0
        end
        norm_convergence = update_P(hf;Δ=Δ,α=0.3)
        
        println("Running HF energy (per moire u.c.): ",Etot)
        println("Running norm convergence: ",norm_convergence)
        push!(iter_energy,Etot)
        push!(iter_err,norm_convergence)
        iter +=1
        
        if iter > 200
            break 
        end
        
    end

    return iter_err, iter_energy
end

function BM_info(hf::HartreeFock)
    """ 
        Add BM kinetic energy term suppressed by parameter α
    """
    tmp_ϵk = reshape(hf.H0,hf.nb*hf.q,hf.nη,hf.ns,hf.nb*hf.q,hf.nη,hf.ns,hf.q,hf.nq,hf.nq)
    tmp_Σz0 = reshape(hf.Σz0,hf.nb*hf.q,hf.nη,hf.ns,hf.nb*hf.q,hf.nη,hf.ns,hf.q,hf.nq,hf.nq)
    hbm = zeros(Float64,hf.nb*hf.q,hf.nq,hf.nq)
    σz = zeros(ComplexF64,hf.nb*hf.q,hf.nb*hf.q,hf.nq,hf.nq)
    
    for iη in 1:2 
        jldopen(hf.metadata[iη]) do file 
            hbm .= file["E"]
            for rk1 in 1:hf.q, iq in 1:(hf.nb*hf.q), is in 1:2
                tmp_ϵk[iq,iη,is,iq,iη,is,rk1,:,:] .= view(hbm,iq,:,:)
            end

            σz .= file["PΣz"]
            for rk1 in 1:hf.q, is in 1:2
                tmp_Σz0[:,iη,is,:,iη,is,rk1,:,:] .=  σz * (3-2iη) 
            end
        end
    end
    # println(maximum(hbm))
    return nothing
end

function add_Hartree(hf::HartreeFock;β::Float64=1.0,V0::Float64=1.0)
    """
        Hartree Contribution suppressed by parameter β
    """
    Lm = sqrt(abs(hf.params.a1)*abs(hf.params.a2))
    tmpΛ = reshape(hf.Λ,2hf.q,hf.nq^2*hf.q,2hf.q,hf.nq^2*hf.q,hf.nη)
    tmpP = reshape(hf.P,2hf.q,hf.nη,hf.ns,2hf.q,hf.nη,hf.ns,hf.nq^2*hf.q)
    tmpH = reshape(hf.H,2hf.q,hf.nη,hf.ns,2hf.q,hf.nη,hf.ns,hf.nq^2*hf.q)
    nk = length(hf.latt.kvec)

    for m in -hf.ng:hf.ng, n in (-hf.ng*hf.q):(hf.ng*hf.q)
        for iη in 1:2
            jldopen(hf.metadata[iη]) do file 
                hf.Λ[:,:,iη] .= file["$(m)_$(n)"]
            end
        end
        trΛ = 0.0 + 0.0im
        if n%hf.q ==0
            for ik in 1:size(tmpΛ,2),is in 1:hf.ns,iη in 1:hf.nη
                trΛ += tr( conj(view(tmpΛ,:,ik,:,ik,iη))*view(tmpP,:,iη,is,:,iη,is,ik) )
            end
        else
            for ik in 1:size(tmpΛ,2),is in 1:hf.ns,iη in 1:hf.nη
                trΛ += tr( conj(view(tmpΛ,:,ik,:,ik,iη))*( view(tmpP,:,iη,is,:,iη,is,ik) + 0.5I) )
            end
        end
        # println(abs(trΛ))
        for ik in 1:size(tmpΛ,2),is in 1:hf.ns,iη in 1:hf.nη
            tmpH[:,iη,is,:,iη,is,ik] .+= (β*V0/nk*V(m*hf.params.g1+n/hf.q*hf.params.g2,Lm)*trΛ) * view(tmpΛ,:,ik,:,ik,iη) 
        end
    end
    return nothing
end

function add_Fock(hf::HartreeFock;β::Float64=1.0,V0::Float64=1.0)
    """
        Fock Contribution 
    """
    Lm = sqrt(abs(hf.params.a1)*abs(hf.params.a2))
    tmpΛ = reshape(hf.Λ,2hf.q,hf.nq^2*hf.q,2hf.q,hf.nq^2*hf.q,hf.nη)
    tmpP = reshape(hf.P,2hf.q,hf.nη,hf.ns,2hf.q,hf.nη,hf.ns,hf.nq^2*hf.q)
    tmpH = reshape(hf.H,2hf.q,hf.nη,hf.ns,2hf.q,hf.nη,hf.ns,hf.nq^2*hf.q)
    nk = length(hf.latt.kvec)
    
    kvec = reshape( reshape(collect(0:(hf.q-1))./hf.q*hf.params.g1,:,1,1) .+ 
                    reshape(hf.latt.k1[1:hf.nq]*hf.params.g1,1,:,1) .+ 
                    reshape(hf.latt.k2[1:hf.nq]*hf.params.g2,1,1,:), : )

    for m in -hf.ng:hf.ng, n in (-hf.ng*hf.q):(hf.ng*hf.q)
        for iη in 1:2
            jldopen(hf.metadata[iη]) do file 
                hf.Λ[:,:,iη] .= file["$(m)_$(n)"]
            end
        end
        for ik in 1:size(tmpΛ,2), iη in 1:hf.nη, is in 1:hf.ns, iηr in 1:hf.nη, isr in 1:hf.ns
            _k = kvec[ik]
            for ip in 1:size(tmpΛ,2)
                _p = kvec[ip]
                tmpH[:,iη,is,:,iηr,isr,ik] .-= (β/nk*V0*V(_p-_k+m*hf.params.g1+n/hf.q*hf.params.g2,Lm)) * 
                        view(tmpΛ,:,ik,:,ip,iη)*transpose(view(tmpP,:,iηr,isr,:,iη,is,ip))*view(tmpΛ,:,ik,:,ip,iηr)'
            end
        end
        # equivalent Fock term: 
        # for ik in 1:size(tmpΛ,2), iη in 1:hf.nη, is in 1:hf.ns, iηr in 1:hf.nη, isr in 1:hf.ns
        #     _k = kvec[ik]
        #     for ip in 1:size(tmpΛ,2)
        #         _p = kvec[ip]
        #         tmpH[:,iη,is,:,iηr,isr,ik] .-= (β/nk*V0*V(_k-_p+m*hf.params.g1+n/hf.q*hf.params.g2,Lm)) * 
        #                 view(tmpΛ,:,ip,:,ik,iη)'*transpose(view(tmpP,:,iηr,isr,:,iη,is,ip))*view(tmpΛ,:,ip,:,ik,iηr)
        #     end
        # end
    end
    return nothing
end


function add_Fock_vectorize(hf::HartreeFock;β::Float64=1.0,V0::Float64=1.0)
    """
        Fock Contribution 
    """
    Lm = sqrt(abs(hf.params.a1)*abs(hf.params.a2))
    tmpΛ = reshape(hf.Λ,2hf.q,hf.nq^2*hf.q,2hf.q,hf.nq^2*hf.q,hf.nη,1)
    tmpP = reshape(hf.P,2hf.q,hf.nη,hf.ns,2hf.q,hf.nη,hf.ns,hf.nq^2*hf.q)
    tmpH = reshape(hf.H,2hf.q,hf.nη,hf.ns,2hf.q,hf.nη,hf.ns,hf.nq^2*hf.q)
    nk = length(hf.latt.kvec)

    kvec = reshape( reshape(collect(0:(hf.q-1))./hf.q*hf.params.g1,:,1,1) .+ 
                    reshape(hf.latt.k1[1:hf.nq]*hf.params.g1,1,:,1) .+ 
                    reshape(hf.latt.k2[1:hf.nq]*hf.params.g2,1,1,:), : )

    P_perm = permutedims(tmpP,(4,1,7,5,6,2,3))
    Λ_perm = zeros(ComplexF64,hf.nb*hf.q,hf.q*hf.nq^2,hf.nb*hf.q,hf.nη,1,hf.q*hf.nq^2)
    Vnum = zeros(Float64,hf.q*hf.nq^2,hf.q*hf.nq^2)

    term1 = zeros(ComplexF64,1,1,hf.q*hf.nq^2,1,1,1,1,1,1,hf.q*hf.nq^2)
    term2 = zeros(ComplexF64,hf.nb*hf.q,1,hf.q*hf.nq^2,hf.nb*hf.q,hf.nη,1,1,1,1,hf.nq^2*hf.q)
    term3 = zeros(ComplexF64,hf.nb*hf.q,hf.nb*hf.q,hf.q*hf.nq^2,1,hf.nη,hf.ns,1,hf.nη,hf.ns,1)
    term4 = zeros(ComplexF64,1,hf.nb*hf.q,hf.q*hf.nq^2,1,1,1,hf.nb*hf.q,hf.nη,1,hf.nq^2*hf.q)
    # vectorize construct matrix: 
    ## H(b,α,p,a,η,s,β,η',s',k), then sum over the first three dimensions
    # size in MB: 256q^6*nq^4*16/1024^2
    # H_to_sum = zeros(ComplexF64,hf.nb*hf.q,hf.nb*hf.q,hf.nq^2*hf.q,
    #                             hf.nb*hf.q,hf.nη,hf.ns,
    #                             hf.nb*hf.q,hf.nη,hf.ns, hf.q*hf.nq^2)
    for m in -hf.ng:hf.ng, n in (-hf.ng*hf.q):(hf.ng*hf.q)
        for iη in 1:2
            jldopen(hf.metadata[iη]) do file 
                hf.Λ[:,:,iη] .= file["$(m)_$(n)"]
            end
        end
        # println(m," ",n," ",norm(hf.Λ))
        Vnum .= (β/nk*V0) .* V.(kvec .-transpose(kvec) .+m*hf.params.g1.+n/hf.q*hf.params.g2,Lm)
        Λ_perm .= permutedims(tmpΛ,(3,4,1,5,6,2))
        # H_to_sum .= reshape(Vnum,(1,1,hf.q*hf.nq^2,1,1,1,1,1,1,hf.q*hf.nq^2)) .* 
        #             reshape(Λ_perm,(hf.nb*hf.q,1,hf.q*hf.nq^2,hf.nb*hf.q,hf.nη,1,1,1,1,hf.nq^2*hf.q)) .* 
        #             reshape(P_perm,(hf.nb*hf.q,hf.nb*hf.q,hf.q*hf.nq^2,1,hf.nη,hf.ns,1,hf.nη,hf.ns,1)) .*
        #             reshape(conj(Λ_perm),(1,hf.nb*hf.q,hf.q*hf.nq^2,1,1,1,hf.nb*hf.q,hf.nη,1,hf.nq^2*hf.q))
        # tmpH .-= reshape(sum( H_to_sum ,dims=(1,2,3)),size(tmpH))
        term1 .= reshape(Vnum,(1,1,hf.q*hf.nq^2,1,1,1,1,1,1,hf.q*hf.nq^2)) 
        term2 .= reshape(Λ_perm,(hf.nb*hf.q,1,hf.q*hf.nq^2,hf.nb*hf.q,hf.nη,1,1,1,1,hf.nq^2*hf.q))
        term3 .= reshape(P_perm,(hf.nb*hf.q,hf.nb*hf.q,hf.q*hf.nq^2,1,hf.nη,hf.ns,1,hf.nη,hf.ns,1))
        term4 .= reshape(conj(Λ_perm),(1,hf.nb*hf.q,hf.q*hf.nq^2,1,1,1,hf.nb*hf.q,hf.nη,1,hf.nq^2*hf.q))
        @inbounds @fastmath for i3 in 1:hf.q*hf.nq^2, i2 in 1:hf.nb*hf.q, i1 in 1:hf.nb*hf.q 
            tmpH .-= view(term1,1,1,i3,:,:,:,:,:,:,:) .* 
                    view(term2,i1,1,i3,:,:,:,:,:,:,:) .*  
                    view(term3,i1,i2,i3,:,:,:,:,:,:,:) .* 
                    view(term4,1,i2,i3,:,:,:,:,:,:,:)
        end
    end
    return nothing
end

function update_P(hf::HartreeFock;Δ::Float64=0.0,α::Float64=0.2)
    """
        Diagonalize Hamiltonian for every k; use ν to keep the lowest N particle states;
        update P 
    """
    νnorm = round(Int,(hf.ν+4)/8 * size(hf.H,1)*size(hf.H,3))  # total number of occupied states 
    vecs = zeros(ComplexF64,size(hf.H,1),size(hf.H,2),size(hf.H,3))
    vals = zeros(Float64,size(hf.H,1),size(hf.H,3))
    for ik in 1:size(hf.H,3)
        check_Hermitian(hf.H[:,:,ik])
        hf.ϵk[:,ik] = eigvals(Hermitian(view(hf.H,:,:,ik)))
        vals[:,ik],vecs[:,:,ik] = eigen(Hermitian(view(hf.H,:,:,ik)-Δ*(conj.(view(hf.P,:,:,ik))+0.5*I)) )
        for iq in 1:size(hf.H,1)
            hf.σzτz[iq,ik] = real(view(vecs,:,iq,ik)'*view(hf.Σz0,:,:,ik)*view(vecs,:,iq,ik))
        end
        # check_Unitary(vecs[:,:,ik])
    end
    # plot_spectra(hf)

    iϵ_sorted = sortperm(vals[:])
    iϵ_occupied = iϵ_sorted[1:νnorm]
    iband_occupied = (iϵ_occupied .-1) .% size(vals,1) .+1
    ik_occupied = (iϵ_occupied .-1) .÷ size(vals,1) .+1

    P_new = zeros(ComplexF64,size(hf.P,1),size(hf.P,2),size(hf.P,3))
    for ik in 1:size(hf.P,3)
        occupied_vecs = vecs[:,iband_occupied[ik_occupied.==ik],ik]
        P_new[:,:,ik] = conj(occupied_vecs)*transpose(occupied_vecs) - 0.5*I
    end

    # α = compute_optimal_λ(hf.H-hf.H0,hf.H0,hf.P,P_new,hf.ν)
    norm_convergence = norm(P_new .- hf.P) ./ norm(P_new)
    hf.P .= α*P_new .+(1-α)*hf.P
    println("ODA value α is: ",α)
    # println("Trace of P is: ", sum([tr(P_new[:,:,ik]) for ik in 1:size(P,3)]) )
    return norm_convergence
end

function compute_optimal_λ(H_HF::Array{ComplexF64,3},H0::Array{ComplexF64,3},P1::Array{ComplexF64,3},P2::Array{ComplexF64,3},ν::Float64)
    # assuming a simple parabolic dependence on λ
    λ = 0.9
    step_λ = 0.02
    E_tot0 = compute_HF_energy(H_HF,H0,(1-λ)*P1+λ*P2,ν)
    while λ>0.1 
        λ -= step_λ
        E_tot = compute_HF_energy(H_HF,H0,(1-λ)*P1+λ*P2,ν)
        if E_tot < E_tot0 
            E_tot0 = E_tot 
        else 
            break 
        end
    end
    return λ
end

function compute_HF_energy(H_HF::Array{ComplexF64,3},H0::Array{ComplexF64,3},P::Array{ComplexF64,3},ν::Float64)
    Etot = 0.0 
    for ik in 1:size(H0,3) 
        Etot += tr(view(H_HF,:,:,ik)*transpose(view(P,:,:,ik)))/2 + tr(view(H0,:,:,ik)*transpose(view(P,:,:,ik)))
    end
    return real(Etot)/(size(H0,3)*size(H0,1))*8 # per moire unit cell
end

function check_Hermitian(H::Matrix{ComplexF64})
    err = norm(H - H')
    if err > 1e-6
        println("Error with Hermitian Hamiltonian")
    end
    return nothing
end


function check_Unitary(A::Matrix{ComplexF64})
    err = norm(A'*A -I)
    if err > 1e-6
        println("Error with Unitarity of Matrix")
    end
    return nothing
end


include("initP_helpers.jl")
# include("plot_helpers.jl")
