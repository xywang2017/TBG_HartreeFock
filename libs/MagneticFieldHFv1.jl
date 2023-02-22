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
    nt::Int
    lk::Int

    # gvec and overlap matrix for a given gvec
    gvec::Matrix{ComplexF64}
    Λ::Array{ComplexF64,4}
    Λs::Array{ComplexF64,4} # prestore all the data of overlap

    # Coulomb unit 
    V0::Float64

    ν::Float64 # filling fraction -4 to 4 
    P::Array{ComplexF64,3}  # one particle density matrix
    ϵk::Matrix{Float64} # eigenvalues of HF renormalized band dispersions nt x ns x nk
    σzτz::Matrix{Float64} # eigenvalues of σzτz for all HF states;

    μ::Float64 # running Hartree-Fock chemical potential 
    Δ::Vector{Float64} # spin-valley-band mixing order parameter s_iη_jn_k
    Δstr::Vector{String} # string

    Σz0::Array{ComplexF64,3} # projected Σz operator in the BM basis
    H0::Array{ComplexF64,3}  # non-interacting part of the Hamiltonian (BM)
    H::Array{ComplexF64,3} # running Hamiltonian for any given k; nt*q x nt*q x nk
    precision::Float64  # iteration stopping point

    savename::String

    HartreeFock() = new()
end

@inline function V(q::ComplexF64,Lm::Float64) ::Float64
    res = 1e-6
    ϵr = 15.0
    return ( abs(q) < res ) ? 0 : 2π/(ϵr*abs(q))*tanh(abs(q)*4Lm/2)
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
        ν::Float64=0.0,ϕ::Rational{Int}=1//10,prefix::String="",_Init::String="CNP",savename::String="placeholder.txt",
        P0::Array{ComplexF64,3}=ones(ComplexF64,1,1,1),H0::Array{ComplexF64,3}=ones(ComplexF64,1,1,1))
    p, q = numerator(ϕ), denominator(ϕ)
    hf.p = p
    hf.q = q
    hf.ν = ν 
    hf.precision = precision
    hf.nb, hf.nη, hf.ns, hf.nt = 2, 2, 2, 8*hf.q # data stored as 2q x nη x ns x nk
    hf.ng = 3
    hf.nq = (q>4) ? 1 : 2
    hf.metadata = [prefix*"_$(p)_$(q)_K_metadata.jld2",
                   prefix*"_$(p)_$(q)_Kprime_metadata.jld2"]
    hf.lk = hf.q*hf.nq^2
    
    hf.savename = savename
    
    hf.params = params
    hf.latt = Lattice() 
    constructLattice(hf.latt,hf.params;lk = hf.nq*hf.q)  #[0,1)x[0,1), so far works for p/q < 1
    hf.V0 = CoulombUnit(hf.params)

    # cannot save the overlap matrix for all of gvecs, so instead only save one gvec at a time
    hf.gvec = reshape(collect(-hf.ng:hf.ng),:,1)*hf.params.g1 .+ reshape(collect(-hf.ng*hf.q:hf.ng*hf.q),1,:)*hf.params.g2 ./hf.q
    hf.Λ = zeros(ComplexF64,hf.nt,hf.lk,hf.nt,hf.lk)

    # BM band structure and projected sublattice info
    hf.H0 = zeros(ComplexF64,hf.nt,hf.nt,hf.lk)
    hf.Σz0 = zeros(ComplexF64,size(hf.H0))
    hf.Λs = zeros(ComplexF64,hf.nb*hf.lk,hf.nb*hf.lk,hf.nη,length(hf.gvec))
    BM_info(hf)

    # order parameters 
    ηs = ["η0","η1","η2","η3"]
    σs = ["s0","s1","s2","s3"]
    hf.Δstr = [σs[k]*ηs[j] for k in 1:4 for j in 1:4]
    hf.Δ = zeros(Float64,size(hf.Δstr))
    # ------------------------------------------------- Begin Hartree Fock Procedure -------------------------------- #
    
    hf.P = zeros(ComplexF64,size(hf.H0))
    hf.H = zeros(ComplexF64,size(hf.H0))
    hf.ϵk = zeros(Float64,size(hf.H0,1),size(hf.H0,3))
    hf.σzτz = zeros(Float64,size(hf.H0,1),size(hf.H0,3))
    
    # --------- Initialization ---------- #
    
    init_P(hf,_Init=_Init,P0=P0,H0=H0)
    
    # --------- Hartree Fock Iterations ---------- #
    norm_convergence = 10.0 
    iter = 1
    iter_err = Float64[]
    iter_energy = Float64[]
    iter_oda = Float64[]
    while norm_convergence > hf.precision
        @time begin 
            println("Iter: ",iter)
            hf.H .= hf.H0 * 1.0
            # add_Hartree(hf;β=1.0)
            # add_Fock(hf;β=1.0)
            add_HartreeFock(hf;β=1.0)
            Etot = compute_HF_energy(hf.H .- hf.H0,hf.H0,hf.P,hf.ν)
            #Δ is a projector to make it closed shell
            if norm_convergence <1e-4
                Δ = 0.0
            else 
                Δ = 0.0
            end
            norm_convergence,λ = update_P(hf;Δ=Δ)
            
            println("Running HF energy (per moire u.c.): ",Etot)
            println("Running norm convergence: ",norm_convergence)
            println("Running ODA paramter λ: ",λ)
        end
        push!(iter_energy,Etot)
        push!(iter_err,norm_convergence)
        push!(iter_oda,λ)
        # if iter ==1 
        #     save("typical_starting_point.jld2","hf",hf,
        #             "iter_err",iter_err,"iter_energy",iter_energy,"iter_oda",iter_oda)
        # end

        if mod(iter,10) ==0 || norm_convergence < hf.precision
            save(hf.savename,"hf",hf,
                    "iter_err",iter_err,"iter_energy",iter_energy,"iter_oda",iter_oda)
        end

        iter +=1
        
        if iter > 300 || λ < 1e-3
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
    tmp_Λs = reshape(hf.Λs,hf.nb*hf.q*hf.lk,hf.nb*hf.q*hf.lk,hf.nη,2hf.ng,:)
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

            for m in -hf.ng:hf.ng, n in (-hf.ng*hf.q):hf.ng*hf.q
                for iη in 1:2
                    jldopen(hf.metadata[iη]) do file 
                        metadata .= file["$(m)_$(n)"]
                        tmpΛs[:,:,iη,m+hf.ng+1,n+hf.ng*hf.q+1] = file["$(m)_$(n)"]
                    end
                end
            end
        end
    end
    # println(maximum(hbm))
    return nothing
end

function add_Hartree(hf::HartreeFock;β::Float64=1.0)
    """
        Hartree Contribution suppressed by parameter β
    """
    Lm = sqrt(abs(hf.params.a1)*abs(hf.params.a2))
    tmpΛ = reshape(hf.Λ,2hf.q,hf.nη,hf.ns,hf.lk,2hf.q,hf.nη,hf.ns,hf.lk)
    metadata = zeros(ComplexF64,2hf.q*hf.lk,2hf.q*hf.lk)
    tmp_metadata = reshape(metadata,2hf.q,hf.lk,2hf.q,hf.lk)
    for ig in eachindex(hf.gvec)
        G = hf.gvec[ig]
        for iη in 1:2
            metadata .= view(hf.Λs,:,:,iη,ig)
            for is in 1:2
                tmpΛ[:,iη,is,:,:,iη,is,:] .= tmp_metadata 
            end
        end
        trPG = 0.0+0.0im
        for ik in 1:size(hf.P,3)
            trPG += tr(view(hf.P,:,:,ik)*conj(view(hf.Λ,:,ik,:,ik)))
        end
        for ik in 1:size(hf.P,3) 
            hf.H[:,:,ik] .+= ( β/hf.latt.nk*hf.V0*V(G,Lm) * trPG) * view(hf.Λ,:,ik,:,ik)
        end
    end
    return nothing
end

function add_Fock(hf::HartreeFock;β::Float64=1.0)
    """
        Fock Contribution 
    """
    Lm = sqrt(abs(hf.params.a1)*abs(hf.params.a2))
    tmpΛ = reshape(hf.Λ,2hf.q,hf.nη,hf.ns,hf.lk,2hf.q,hf.nη,hf.ns,hf.lk)
    kvec = reshape( reshape(collect(0:(hf.q-1))./hf.q*hf.params.g1,:,1,1) .+ 
                    reshape(hf.latt.k1[1:hf.nq]*hf.params.g1,1,:,1) .+ 
                    reshape(hf.latt.k2[1:hf.nq]*hf.params.g2,1,1,:), : )
    tmp_Fock = zeros(ComplexF64,hf.nt,hf.nt)
    
    metadata = zeros(ComplexF64,2hf.q*hf.lk,2hf.q*hf.lk)
    tmp_metadata = reshape(metadata,2hf.q,hf.lk,2hf.q,hf.lk)

    for ig in eachindex(hf.gvec)
        G = hf.gvec[ig]
        for iη in 1:2
            metadata .= view(hf.Λs,:,:,iη,ig)
            for is in 1:2
                tmpΛ[:,iη,is,:,:,iη,is,:] .= tmp_metadata 
            end
        end
        for ik in 1:size(hf.P,3)
            tmp_Fock .= 0.0 + 0.0im
            for ip in 1:size(hf.P,3)
                tmp_Fock .+= ( β*hf.V0*V(kvec[ip]-kvec[ik]+G,Lm) /hf.latt.nk) * 
                            ( view(hf.Λ,:,ik,:,ip)*transpose(view(hf.P,:,:,ip))*view(hf.Λ,:,ik,:,ip)' )
            end
            hf.H[:,:,ik] .-= tmp_Fock
        end
    end
    return nothing
end

function add_HartreeFock(hf::HartreeFock;β::Float64=1.0)
    Lm = sqrt(abs(hf.params.a1)*abs(hf.params.a2))
    tmpΛ = reshape(hf.Λ,2hf.q,hf.nη,hf.ns,hf.lk,2hf.q,hf.nη,hf.ns,hf.lk)
    kvec = reshape( reshape(collect(0:(hf.q-1))./hf.q*hf.params.g1,:,1,1) .+ 
                    reshape(hf.latt.k1[1:hf.nq]*hf.params.g1,1,:,1) .+ 
                    reshape(hf.latt.k2[1:hf.nq]*hf.params.g2,1,1,:), : )
    tmp_Fock = zeros(ComplexF64,hf.nt,hf.nt)
    
    metadata = zeros(ComplexF64,2hf.q*hf.lk,2hf.q*hf.lk)
    tmp_metadata = reshape(metadata,2hf.q,hf.lk,2hf.q,hf.lk)

    for ig in eachindex(hf.gvec)
        G = hf.gvec[ig]
        for iη in 1:2
            metadata .= view(hf.Λs,:,:,iη,ig)
            for is in 1:2
                tmpΛ[:,iη,is,:,:,iη,is,:] .= tmp_metadata 
            end
        end
        # --------------------------------------- Hartree ------------------------------- #
        trPG = 0.0+0.0im
        for ik in 1:size(hf.P,3)
            trPG += tr(view(hf.P,:,:,ik)*conj(view(hf.Λ,:,ik,:,ik)))
        end
        for ik in 1:size(hf.P,3) 
            hf.H[:,:,ik] .+= ( β/hf.latt.nk*hf.V0*V(G,Lm) * trPG) * view(hf.Λ,:,ik,:,ik)
        end
        # --------------------------------------- Fock ------------------------------- #
        for ik in 1:size(hf.P,3)
            tmp_Fock .= 0.0 + 0.0im
            for ip in 1:size(hf.P,3)
                tmp_Fock .+= ( β*hf.V0*V(kvec[ip]-kvec[ik]+G,Lm) /hf.latt.nk) * 
                            ( view(hf.Λ,:,ik,:,ip)*transpose(view(hf.P,:,:,ip))*view(hf.Λ,:,ik,:,ip)' )
            end
            hf.H[:,:,ik] .-= tmp_Fock
        end
    end
    return nothing
end

function update_P(hf::HartreeFock;Δ::Float64=0.0)
    """
        Diagonalize Hamiltonian for every k; use ν to keep the lowest N particle states;
        update P 
    """
    νnorm = round(Int,(hf.ν+4)/8 * size(hf.H,1)*size(hf.H,3))  # total number of occupied states 
    vecs = zeros(ComplexF64,size(hf.H,1),size(hf.H,2),size(hf.H,3))
    vals = zeros(Float64,size(hf.H,1),size(hf.H,3))
    for ik in 1:size(hf.H,3)
        # check_Hermitian(hf.H[:,:,ik])
        hf.ϵk[:,ik] = eigvals(Hermitian(view(hf.H,:,:,ik)))
        vals[:,ik],vecs[:,:,ik] = eigen(Hermitian(view(hf.H,:,:,ik)-Δ*(conj.(view(hf.P,:,:,ik))+0.5*I)) )
        for iq in 1:size(hf.H,1)
            hf.σzτz[iq,ik] = real(view(vecs,:,iq,ik)'*view(hf.Σz0,:,:,ik)*view(vecs,:,iq,ik))
        end
        # check_Unitary(vecs[:,:,ik])
    end
    # plot_spectra(hf)

    hf.μ = find_chemicalpotential(hf.ϵk[:],νnorm/(size(hf.H,1)*size(hf.H,3)))
    hf.Δ .= calculate_valley_spin_band_order_parameters(hf)

    iϵ_sorted = sortperm(vals[:])
    iϵ_occupied = iϵ_sorted[1:νnorm]
    iband_occupied = (iϵ_occupied .-1) .% size(vals,1) .+1
    ik_occupied = (iϵ_occupied .-1) .÷ size(vals,1) .+1

    P_new = zeros(ComplexF64,size(hf.P,1),size(hf.P,2),size(hf.P,3))
    for ik in 1:size(hf.P,3)
        occupied_vecs = vecs[:,iband_occupied[ik_occupied.==ik],ik]
        P_new[:,:,ik] = conj(occupied_vecs)*transpose(occupied_vecs) - 0.5*I
    end

    λ = oda_parametrization(hf,P_new .- hf.P;β=1.0)
    # λ = 1.0 # often times oda_parameterization returns λ = 1.0, therefore not necessary
    norm_convergence = calculate_norm_convergence(λ*P_new + (1-λ)*hf.P,hf.P)
    hf.P .= λ*P_new + (1-λ)*hf.P
    return norm_convergence,λ
end


function calculate_norm_convergence(P2::Array{ComplexF64,3},P1::Array{ComplexF64,3})
    # vals1 = zeros(Float64,size(P1,1),size(P1,3))
    # vals2 = zeros(Float64,size(P1,1),size(P1,3))
    # for ik in 1:size(P1,3)
    #     vals1[:,ik] .= eigvals(Hermitian(view(P1,:,:,ik)))
    #     vals2[:,ik] .= eigvals(Hermitian(view(P2,:,:,ik))) 
    # end
    # return norm(vals2 .-vals1) / norm(vals2)
    return norm(P1 .- P2) ./ norm(P2)
end

function compute_HF_energy(H_HF::Array{ComplexF64,3},H0::Array{ComplexF64,3},P::Array{ComplexF64,3},ν::Float64)
    Etot = 0.0 
    for ik in 1:size(H0,3) 
        Etot += tr(view(H_HF,:,:,ik)*transpose(view(P,:,:,ik)))/2 + tr(view(H0,:,:,ik)*transpose(view(P,:,:,ik)))
    end
    return real(Etot)/(size(H0,3)*size(H0,1))*8 # per moire unit cell
end

function oda_parametrization(hf::HartreeFock,δP::Array{ComplexF64,3};β::Float64=1.0)
    # compute coefficients b λ + a λ^2/2

    Lm = sqrt(abs(hf.params.a1)*abs(hf.params.a2))
    tmpΛ = reshape(hf.Λ,2hf.q,hf.nη,hf.ns,hf.lk,2hf.q,hf.nη,hf.ns,hf.lk)
    kvec = reshape( reshape(collect(0:(hf.q-1))./hf.q*hf.params.g1,:,1,1) .+ 
                    reshape(hf.latt.k1[1:hf.nq]*hf.params.g1,1,:,1) .+ 
                    reshape(hf.latt.k2[1:hf.nq]*hf.params.g2,1,1,:), : )
    tmp_Fock = zeros(ComplexF64,hf.nt,hf.nt)

    metadata = zeros(ComplexF64,2hf.q*hf.lk,2hf.q*hf.lk)
    tmp_metadata = reshape(metadata,2hf.q,hf.lk,2hf.q,hf.lk)

    # change of Hartree-Fock due to a small δP
    δH = zeros(ComplexF64,size(δP))
    for ig in eachindex(hf.gvec)
        G = hf.gvec[ig]
        for iη in 1:2
            metadata .= view(hf.Λs,:,:,iη,ig)
            for is in 1:2
                tmpΛ[:,iη,is,:,:,iη,is,:] .= tmp_metadata 
            end
        end
        trPG = 0.0+0.0im
        for ik in 1:size(δH,3) 
            trPG += tr(view(δP,:,:,ik)*conj(view(hf.Λ,:,ik,:,ik)))
        end
        for ik in 1:size(δH,3) 
            δH[:,:,ik] .+= ( β/hf.latt.nk*hf.V0*V(G,Lm) * trPG) * view(hf.Λ,:,ik,:,ik)
        end
        for ik in 1:size(δH,3) 
            tmp_Fock .= 0.0 + 0.0im
            for ip in 1:size(δH,3) 
                tmp_Fock .+= ( β*hf.V0*V(kvec[ip]-kvec[ik]+G,Lm) /hf.latt.nk) * 
                            ( view(hf.Λ,:,ik,:,ip)*transpose(view(δP,:,:,ip))*view(hf.Λ,:,ik,:,ip)' )
            end
            δH[:,:,ik] .-= tmp_Fock
        end
    end

    # compute coefficients with δP 
    a, b = 0.0+0.0im , 0.0 +0.0im
    for ik in 1:size(δH,3) 
        b += tr(transpose(view(δP,:,:,ik))*view(hf.H0,:,:,ik)) + 
             tr(transpose(view(δP,:,:,ik))*view(hf.H .- hf.H0,:,:,ik))/2 +
             tr(transpose(view(hf.P,:,:,ik))*view(δH,:,:,ik))/2
        a += tr(transpose(view(δP,:,:,ik))*view(δH,:,:,ik))
    end
    a = real(a)/size(δP,3)
    b = real(b)/size(δP,3)

    λ, λ0 = 0.0 , -b/a
    if a>0 # convex and increasing with large λ 
        if λ0 <=0 
            λ = 0.1  # give it some kick..
        elseif λ0 <1 
            λ = λ0 
        else
            λ = 1.0 
        end
    else
        if λ0 <=0.5 
            λ = 1.0 
        else 
            λ = 0.1 # give it some kick..
        end
    end
    return λ
end


function calculate_valley_spin_band_order_parameters(hf::HartreeFock)
    s0 = ComplexF64[1 0;0 1]
    s1 = ComplexF64[0 1;1 0]
    s2 = ComplexF64[0 -1im;1im 0]
    s3 = ComplexF64[1 0;0 -1]
    pauli_matrices = [s0,s1,s2,s3]
    I_band = Array{Float64,2}(I,hf.nb*hf.q,hf.nb*hf.q)
    order_parameters = Float64[]
    Δ = zeros(Float64,size(hf.ϵk))
    for k in 1:4
        for j in 1:4 
            for ik in 1:size(hf.ϵk,2)
                F = eigen(Hermitian(view(hf.H,:,:,ik)))
                Δ[:,ik] = real(diag(F.vectors'*kron(pauli_matrices[k],kron(pauli_matrices[j],I_band))*F.vectors))
            end
            push!(order_parameters,sum(Δ[:][hf.ϵk[:].<= hf.μ])/length(hf.ϵk)*8)
        end
    end
    return order_parameters
end

function find_chemicalpotential(energies::Vector{Float64},ν::Float64)
    E = sort(energies)
    νs = collect(eachindex(E)) ./ length(E)
    i = 1
    while ν > νs[i]
        i = i+1 
    end   
    if i < length(E)
        μ = (E[i+1]+E[i])/2
    else
        μ = E[i]
    end
    return μ
end

include("initP_helpers.jl")
# include("plot_helpers.jl")
