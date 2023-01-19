using Arpack
include("Parameters_mod.jl")
include("Lattice_mod.jl")
# --------------------------------------------------------------------------------------------------------------- #
mutable struct HBM
    flag_valley::String # if "K" only calculate K valey; if "Both" construct opposite valley by symmetry 
    _σrotation::Bool # flag on whether to add σ matrix rotation 

    params::Params 
    latt::Lattice
    
    nη::Int # number of valleys 
    ns::Int # number of spins
    nb::Int # number of bands 
    nfl::Int # total number of flavors

    hbm::Vector{Float64} # Hbm energies nfl x nk
    Λkp::Array{ComplexF64,3} # αk x βp x G 
    Gs::Vector{ComplexF64} # Gs for Λkp

    HBM() = new()
end


function initHBM(blk::HBM,latt::Lattice,params::Params;lg::Int=9,_σrotation::Bool=true)
    blk._σrotation = _σrotation
    blk.params = params 
    blk.latt = latt
    blk.flag_valley = "Both"
    blk.nη = 2
    blk.ns = 2
    blk.nb = 2 # flat bands 
    blk.nfl = blk.nη*blk.nb # treat spin separately 

    s0 = Float64[1 0; 0 1]
    s1 = Float64[0 1; 1 0]
    is2 = Float64[0 1; -1 0]

    @assert (lg-1)%2 == 0   # we deal with lg being odd, such that -g and g are related easily
    
    gvec = zeros(ComplexF64,lg^2)
    for ig in 1:lg^2
        ig1, ig2 = (ig-1)%lg+1 -(lg+1)÷2 , (ig-1)÷lg+1 - (lg+1)÷2
        gvec[ig] = params.g1 * ig1 + params.g2 * ig2
    end

    # this gives i mu_y I operation in the Bloch basis
    Ig = reverse(Array{Float64}(I,lg^2,lg^2),dims=1)
    Ph = -kron(Ig,kron(is2,s0))

    # this gives C2T eigenstates
    Ig = Array{Float64}(I,lg^2,lg^2)
    C2T = kron(Ig,kron(s0,s1)) # × conj(...)

    nlocal = 4 # layer x sublattice

    T12 = zeros(ComplexF64,nlocal*lg^2,nlocal*lg^2)
    generate_T12(T12,lg,params)

    # temporary container H for each k
    H = zeros(ComplexF64,nlocal*lg^2,nlocal*lg^2)
    # temporary basis vectors and eigenenergies for constructing HBM
    Uk = zeros(ComplexF64,nlocal*lg^2,blk.nb,latt.nk)
    Hk =zeros(Float64,blk.nb,latt.nk)
    for ik in 1:latt.nk
        kval = real(latt.kvec[ik])*params.g1 + imag(latt.kvec[ik])*params.g2
        ComputeH(H,kval,gvec,lg,params,_σrotation)
        H .= H + T12 - params.μ*I
        # Find the smallest eigenvalue and eigenvectors close to zero
        vals, vecs = eigs(Hermitian(H),nev=blk.nb,which=:SM)
        # C2T is broken for hBN alignment
        vecs = vecs + C2T*conj(vecs)
        for i in 1:blk.nb
            tmp = view(vecs,:,i)
            normalize!(tmp)
        end

        if (norm(imag(vals))<1e-6)
            perm = sortperm(real(vals[:]))
            Uk[:,:,ik] = view(vecs,:,perm)
            Hk[:,ik] = real(vals[perm])
        else
            print("Error with Hermiticity of Hamiltonian!\n")
        end
    end

    if isequal(blk.flag_valley,"Both")
        blk.hbm = zeros(Float64,blk.nfl*latt.nk)
        tmp_hbm = reshape(blk.hbm,blk.nη,blk.nb,latt.nk)
        tmp_hbm[1,:,:] = Hk
        if latt.flag_inv == true
            tmp_hbm[2,:,:] = Hk[:,latt.nk:(-1):1]
        else
            println("Lattice is not inversion symmetric! Function initBM")
        end
    elseif isequal(blk.flag_valley,"K")
        blk.hbm = zeros(Float64,blk.nb*latt.nk)
        tmp_hbm = reshape(blk.hbm,blk.nb,latt.nk)
        tmp_hbm[:,:] = Hk
    else 
        println("Wrong flag_valley value")
    end

    if true
        # construct Λ matrices 
        lG = 7 
        Gs = collect((-(lG-1)÷2):((lG-1)÷2))
        blk.Gs = (reshape(Gs,:,1)*params.g1 .+ reshape(Gs,1,:)*params.g2)[:]
        λkp = zeros(ComplexF64,blk.nb*latt.nk,blk.nb*latt.nk)
        λkpKprime = zeros(ComplexF64,blk.nb*latt.nk,blk.nb*latt.nk)
        ur = zeros(ComplexF64,nlocal*lg^2,blk.nb*latt.nk)
        ul = reshape(Uk,:,blk.nb*latt.nk)'
        if isequal(blk.flag_valley,"Both")
            blk.Λkp = zeros(ComplexF64,blk.nfl*latt.nk,blk.nfl*latt.nk,lG^2)
            tmp_Λkp = reshape(blk.Λkp,blk.nη,blk.nb*latt.nk,blk.nη,blk.nb*latt.nk,lG^2)
            tmp_Uk = reshape(Uk,nlocal,lg,lg,blk.nb*latt.nk)
            for ig in 1:size(blk.Λkp,3)
                m,n = Gs[(ig-1)%lG+1],Gs[(ig-1)÷lG+1]
                ur .= reshape(circshift(tmp_Uk,(0,-m,-n,0)),nlocal*lg^2,blk.nb*latt.nk)
                λkp .= ul * ur
                tmp_Λkp[1,:,1,:,ig] = λkp
                λkpKprime .= reshape(reshape(λkp,blk.nb,latt.nk,blk.nb,latt.nk)[:,latt.nk:(-1):1,:,latt.nk:(-1):1],blk.nb*latt.nk,blk.nb*latt.nk)
                tmp_Λkp[2,:,2,:,lG^2-ig+1] = λkpKprime
            end
        elseif isequal(blk.flag_valley,"K")
            blk.Λkp = zeros(ComplexF64,blk.nb*latt.nk,blk.nb*latt.nk,lG^2)
            tmp_Uk = reshape(Uk,nlocal,lg,lg,blk.nb*latt.nk)
            for ig in 1:size(blk.Λkp,3)
                m,n = Gs[(ig-1)%lG+1],Gs[(ig-1)÷lG+1]
                ur .= reshape(circshift(tmp_Uk,(0,-m,-n,0)),nlocal*lg^2,blk.nb*latt.nk)
                λkp .= ul * ur
                blk.Λkp[:,:,ig] = λkp
            end
        else 
            println("Wrong flag_valley value")
        end
    end

    return nothing
end

@inline function dirac(k::ComplexF64,θ0::Float64) ::Matrix{ComplexF64}
    return  abs(k)*[0 exp(-1im*(angle(k)-θ0));exp(1im*(angle(k)-θ0)) 0]
end

function generate_T12(T12::Matrix{ComplexF64},lg::Int,params::Params)
    # p.b.c. is used 
    idg = reshape(collect(1:lg^2),lg,lg)
    # per Oskar & Jian choice of g1 and g2
    idg_nn1 = circshift(idg,(0,1))  # T1 * (|t><b|)
    idg_nn2 = circshift(idg,(1,1))  # T2 * (|t><b|)
    idg_nn12 = circshift(idg,(0,0))  # T0 * (|t><b|)

    tmp = reshape(T12,4,lg^2,4,lg^2)
    for ig in eachindex(idg)
        tmp[3:4,idg[ig],1:2,idg_nn1[ig]] = params.T2
        tmp[1:2,idg_nn1[ig],3:4,idg[ig]] = params.T2'
        tmp[3:4,idg[ig],1:2,idg_nn2[ig]] = params.T1
        tmp[1:2,idg_nn2[ig],3:4,idg[ig]] = params.T1'
        tmp[3:4,idg[ig],1:2,idg_nn12[ig]] = params.T0
        tmp[1:2,idg_nn12[ig],3:4,idg[ig]] = params.T0'
    end
    return nothing
end

function ComputeH(H::Matrix{ComplexF64},k::ComplexF64,gvec::Vector{ComplexF64},lg::Int,params::Params,_σrotation::Bool;ig::Vector{Int}=[0;0])
    """
        Dirac Hamiltonian in the Bloch band basis; Note here k only takes values within first mBZ 
    """
    H .= 0.0 + 0.0im
    σz = ComplexF64[1 0 ; 0 -1]
    σ0 = ComplexF64[1 0; 0 1]
    itr = reshape(collect(1:lg^2),lg,lg)
    idg = view(circshift(itr,(ig[1],ig[2])),:)
    R = params.dθ/2 * Float64[0 -1;1 0]
    ∇u = (params.S[1,1] + params.S[2,2])/2
    # dispersive part
    for ig in 1:lg^2
        qc = gvec[ig]
        kb = k - params.Kb + qc
        kt = k - params.Kt + qc
        if (_σrotation==true)
            k1 = (I + R - params.S*params.α)*[real(kb);imag(kb)]
            k2 = (I - R + params.S*(1-params.α))*[real(kt);imag(kt)]
        else
            k1 = [real(kb);imag(kb)]
            k2 = [real(kt);imag(kt)]
        end
        H[(4idg[ig]-3):(4idg[ig]-2),(4idg[ig]-3):(4idg[ig]-2)] = params.vf*dirac(k1[1]+1im*k1[2],0.0) .- (params.Da * ∇u)*σ0
        H[(4idg[ig]-1):(4idg[ig]),(4idg[ig]-1):(4idg[ig])] = params.vf*dirac(k2[1]+1im*k2[2],0.0) .+ (params.Da * ∇u)*σ0
    end
    return nothing
end
