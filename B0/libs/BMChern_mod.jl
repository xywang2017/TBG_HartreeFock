include("ParametersChern_mod.jl")
include("Lattice_mod.jl")
include("helpers.jl")
# --------------------------------------------------------------------------------------------------------------- #
mutable struct HBM
    _σrotation::Bool # flag on whether to add σ matrix rotation 

    params::Params 
    latt::Lattice
    
    nlocal::Int  # layer x sublattice
    nη::Int # number of valleys 
    ns::Int # number of spins
    nb::Int # number of bands 
    nt::Int # total number of flavors

    H::Array{ComplexF64,4}
    Σz::Array{ComplexF64,3}  #sublattice operator, will save in data file as σzτz
    Uk::Array{ComplexF64,4}  # vectors in chern basis 
    spectrum::Array{Float64,3} # Hbm energies nb x nη x nk
    chern_ham::Array{ComplexF64,4} # BM in chern basis: nb x nb x nη x nk 
    lg::Int 
    gvec::Vector{ComplexF64}


    Λ::Array{ComplexF64,2} # αk x βp 
    Gs::Vector{ComplexF64} # Gs for Λkp

    fname::String # file store name

    HBM() = new()
end

function initHBM(A::HBM,latt::Lattice,params::Params;lg::Int=9,
        _σrotation::Bool=true,_calculate_overlap::Bool=true,fname::String="holder.jld2")
    A._σrotation = _σrotation
    A.params = params 
    A.latt = latt
    A.nη,A.ns,A.nb, A.nlocal = 2, 2, 2, 4
    A.nt = A.nη*A.ns*A.nb 
    A.fname = fname

    @assert (lg-1)%2 == 0   # we deal with lg being odd, such that -g and g are related easily
    A.lg = lg
    A.gvec = ( reshape(collect((-lg÷2):(lg÷2)),:,1)*params.g1 .+ reshape(collect((-lg÷2):(lg÷2)),1,:)*params.g2 )[:]

    constructHBM(A)
    computeSpectrum(A)
    enforceSymmetry(A)
    calculateChern(A)

    jldopen(A.fname, "w") do file
        hbm = zeros(ComplexF64,A.nt,A.nt,A.latt.nk)
        tmp_hbm = reshape(hbm,A.ns,A.nη,A.nb,A.ns,A.nη,A.nb,A.latt.nk)
        for is in 1:A.ns, iη in 1:A.nη
            tmp_hbm[is,iη,:,is,iη,:,:] = A.chern_ham[:,:,iη,:] 
        end
        file["E"] = hbm
        file["nη"] = A.nη
        file["nb"] = A.nb
        file["ns"] = A.ns
        file["Σz"] = A.Σz
    end

    if _calculate_overlap 
        lG = 7
        Gs = collect((-(lG-1)÷2):((lG-1)÷2))
        A.Gs = (reshape(Gs,:,1)*A.params.g1 .+ reshape(Gs,1,:)*A.params.g2)[:]
        A.Λ = zeros(ComplexF64,A.nt*A.latt.nk,A.nt*A.latt.nk)
        jldopen(A.fname,"a") do file
            file["Gs"] = A.Gs
            file["lG"] = lG
            for ig in eachindex(A.Gs)
                m,n = Gs[(ig-1)%lG+1],Gs[(ig-1)÷lG+1]
                println("m: ",m," n:",n)
                calculate_overlap(A,m,n)
                file["$(m)_$(n)"] = A.Λ
            end
        end
    end
    
    return nothing 
end

function constructHBM(A::HBM)
    A.H = zeros(ComplexF64,A.nlocal*A.lg^2,A.nlocal*A.lg^2,A.nη,A.latt.nk)
    T12 = zeros(ComplexF64,A.nlocal*A.lg^2,A.nlocal*A.lg^2)
    H0 = zeros(ComplexF64,size(A.H,1),size(A.H,2))
    for iη in 1:A.nη
        ζ = 3-2iη
        generate_T12(T12,ζ,A)
        for ik in 1:A.latt.nk
            kval = A.latt.kvec[ik]
            constructDiagonals(H0,kval,ζ,A)
            A.H[:,:,iη,ik] = H0 + T12 - params.μ*I
        end
    end
    return nothing 
end 

function computeSpectrum(A::HBM)
    A.spectrum = zeros(Float64,A.nb,A.nη,A.latt.nk)
    A.Uk = zeros(ComplexF64,A.nlocal*A.lg^2,A.nb,A.nη,A.latt.nk)
    
    for ik in 1:A.latt.nk, iη in 1:A.nη
        # check_Hermitian(A.H[:,:,iη,ik])
        H0 = view(A.H,:,:,iη,ik)
        idx_mid = size(H0,1)÷2
        A.spectrum[:,iη,ik], A.Uk[:,:,iη,ik] = eigen(Hermitian(H0),(idx_mid-A.nb÷2+1):(idx_mid+A.nb÷2))
        
    end

    return nothing
end

function enforceSymmetry(A::HBM)
    s0 = Float64[1 0; 0 1]
    s1 = Float64[0 1; 1 0]
    is2 = Float64[0 1; -1 0]
    s3 = Float64[1 0;0 -1]
    # this gives C2T eigenstates
    Ig = Array{Float64}(I,A.lg^2,A.lg^2)
    C2T = kron(Ig,kron(s0,s1)) # × conj(...)

    # if C2T is preserved then xxx. note that it is broken for hBN alignment
    for ik in 1:size(A.Uk,4), iη in 1:A.nη
        A.Uk[:,:,iη,ik] = (view(A.Uk,:,:,iη,ik) .+ C2T*conj.(view(A.Uk,:,:,iη,ik)))
        for j in 1:size(A.Uk,2)
            normalize!(view(A.Uk,:,j,iη,ik))
        end
    end

    # this gives Ph i mu_y I operation in the Bloch basis
    Ig = reverse(Array{Float64}(I,A.lg^2,A.lg^2),dims=1)
    Ph = -kron(Ig,kron(is2,s0))
    for iη in 1:A.nη, ik in 1:A.latt.nk
        # use lower BM band to constrain upper band 
        A.Uk[:,2,iη,ik] = Ph*view(A.Uk,:,1,iη,A.latt.nk+1-ik)
        # val = abs( abs(view(tmpU,:,2,iη,ik)'*view(A.Uk,:,2,iη,ik)) -1)
        # if val > 1e-6
        #     println("error with applying PH: ",val," ")
        # end
    end

    # C2 σ1I 
    # tmpU = deepcopy(A.Uk)
    Ig = reverse(Array{Float64}(I,A.lg^2,A.lg^2),dims=1)
    C2 = kron(Ig,kron(s0,s1))
    for ik in 1:A.latt.nk
        A.Uk[:,:,2,ik] = C2*view(A.Uk,:,:,1,A.latt.nk+1-ik)
        # val = norm( abs.(view(tmpU,:,:,2,ik)'*view(A.Uk,:,:,2,ik)) -I)
        # if val > 1e-6
        #     println("error with applying C2: ",val," ")
        # end
    end

    # going to Chern basis 
    UChern = deepcopy(A.Uk)
    A.chern_ham = zeros(ComplexF64,A.nb,A.nb,A.nη,A.latt.nk)
    vec = zeros(ComplexF64,size(A.Uk,1))
    mat = ComplexF64[1 1;1im -1im]/sqrt(2)
    σz = kron(Array{Float64}(I,A.lg^2,A.lg^2),kron(s0,s3))
    for iη in 1:A.nη, ik in 1:A.latt.nk 
        vec .= (view(A.Uk,:,1,iη,ik) .+ 1im * view(A.Uk,:,2,iη,ik))./sqrt(2)
        if real(vec'*σz*vec)*(3-2iη) < 0
            A.Uk[:,2,iη,ik] = view(A.Uk,:,2,iη,ik) * (-1)
        end
        UChern[:,1,iη,ik] = (view(A.Uk,:,1,iη,ik) + 1im * view(A.Uk,:,2,iη,ik))./sqrt(2)
        UChern[:,2,iη,ik] = (view(A.Uk,:,1,iη,ik) - 1im * view(A.Uk,:,2,iη,ik))./sqrt(2)
        for j in 1:A.nb 
            normalize!(view(UChern,:,j,iη,ik))
        end
        A.chern_ham[:,:,iη,ik] = mat' * Diagonal(A.spectrum[:,iη,ik]) * mat
    end

    # smooth gauge for chern basis along g1 direction
    tmpUChern = reshape(UChern,:,2,2,A.latt.lk,A.latt.lk) 
    for iη in 1:A.nη
        for ik2 in 1:A.latt.lk, ik1 in 2:A.latt.lk
            if real(view(tmpUChern,:,1,iη,ik1-1,ik2)'*view(tmpUChern,:,1,iη,ik1,ik2)) <0 
                tmpUChern[:,:,iη,ik1,ik2] .*= -1 
            end
        end
        for ik2 in 2:A.latt.lk 
            if real(view(tmpUChern,:,1,iη,1,ik2-1)'*view(tmpUChern,:,1,iη,1,ik2)) <0 
                tmpUChern[:,:,iη,:,ik2] .*= -1 
            end
        end  
    end
    A.Uk .= UChern
    # test energies: 
    # for iη in 1:A.nη, ik in 1:A.latt.nk 
    #     vals = eigvals(A.chern_ham[:,:,iη,ik])
    #     norm_diff = norm(vals - A.spectrum[:,iη,ik])
    #     if norm_diff > 1e-6
    #         println("error with Chern hamiltonian: ",norm_diff," ")
    #     end
    # end
    return nothing
end

function calculateChern(A::HBM)
    A.Σz = zeros(ComplexF64,A.ns*A.nη*A.nb,A.ns*A.nη*A.nb,A.latt.nk)
    tmpΣz = reshape(A.Σz,A.ns,A.nη,A.nb,A.ns,A.nη,A.nb,A.latt.nk)

    s0 = ComplexF64[1 0; 0 1]
    sz = ComplexF64[1 0; 0 -1]
    σz = kron(Array{ComplexF64,2}(I,A.lg^2,A.lg^2),kron(s0,sz))

    mat2x2 = zeros(ComplexF64,A.nb,A.nb)
    for ik in 1:A.latt.nk, iη in 1:A.nη
        mat2x2 .= view(A.Uk,:,:,iη,ik)' * σz * view(A.Uk,:,:,iη,ik) * (3-2iη) 
        for is in 1:A.ns
            tmpΣz[is,iη,:,is,iη,:,ik] = mat2x2
        end
    end
    return nothing
end

@inline function dirac(k::ComplexF64,ζ::Int,θ0::Float64) ::Matrix{ComplexF64}
    #ζ = 1 (K) and -1 (K')
    return  ζ*abs(k)*[0 exp(-1im*ζ*(angle(k)-θ0));exp(1im*ζ*(angle(k)-θ0)) 0]
end

function generate_T12(T12::Matrix{ComplexF64},ζ::Int,A::HBM)
    T12 .= 0.0 + 0.0im
    # p.b.c. is used 
    idg = reshape(collect(1:A.lg^2),A.lg,A.lg)
    # # per Oskar & Jian choice of g1 and g2
    # idg_nn1 = circshift(idg,(0,ζ))  # T1 * (|t><b|)
    # idg_nn2 = circshift(idg,(ζ,ζ))  # T2 * (|t><b|)
    # idg_nn12 = circshift(idg,(0,0))  # T0 * (|t><b|)
    # chosing Kt = - Kb, implicitly inversion symmetric
    idg_nn1 = circshift(idg,(-ζ,ζ))  # T1 * (|t><b|)
    idg_nn2 = circshift(idg,(0,ζ))  # T2 * (|t><b|)
    idg_nn12 = circshift(idg,(-ζ,0))  # T0 * (|t><b|)

    tmp = reshape(T12,A.nlocal,A.lg^2,A.nlocal,A.lg^2)
    if ζ ==1 
        T0, T1, T2 = A.params.T0, A.params.T1, A.params.T2 
    elseif ζ == -1
        T0, T1, T2 = A.params.T0, A.params.T2, A.params.T1
    else 
        println("T12: Wrong value of ζ (valley) flag!")
    end

    for ig in eachindex(idg)
        tmp[3:4,idg[ig],1:2,idg_nn1[ig]] = T2
        tmp[1:2,idg_nn1[ig],3:4,idg[ig]] = T2
        tmp[3:4,idg[ig],1:2,idg_nn2[ig]] = T1
        tmp[1:2,idg_nn2[ig],3:4,idg[ig]] = T1
        tmp[3:4,idg[ig],1:2,idg_nn12[ig]] = T0
        tmp[1:2,idg_nn12[ig],3:4,idg[ig]] = T0
    end
    
    return nothing
end

function constructDiagonals(H::Matrix{ComplexF64},k::ComplexF64,ζ::Int,A::HBM)
    """
        Dirac Hamiltonian in the Bloch band basis; Note here k only takes values within first mBZ 
        ζ = 1 (K), -1 (K')
    """
    H .= 0.0 + 0.0im
    σz = ComplexF64[1 0 ; 0 -1]
    σ0 = ComplexF64[1 0; 0 1]
    R = -A.params.dθ/2 * Float64[0 -1;1 0]
    ∇u = (A.params.S[1,1] + A.params.S[2,2])/2
    # dispersive part
    for ig in 1:A.lg^2
        qc = A.gvec[ig]
        kb = k - ζ*A.params.Kb + qc
        kt = k - ζ*A.params.Kt + qc
        if (A._σrotation==true)
            k1 = (I + R - A.params.S*A.params.α)*[real(kb);imag(kb)]
            k2 = (I - R + A.params.S*(1-A.params.α))*[real(kt);imag(kt)]
        else
            k1 = [real(kb);imag(kb)]
            k2 = [real(kt);imag(kt)]
        end
        H[(4ig-3):(4ig-2),(4ig-3):(4ig-2)] = A.params.vf*dirac(k1[1]+1im*k1[2],ζ,0.0) .- (A.params.Da * ∇u)*σ0
        H[(4ig-1):(4ig),(4ig-1):(4ig)] = A.params.vf*dirac(k2[1]+1im*k2[2],ζ,0.0) .+ (A.params.Da * ∇u)*σ0
    end
    return nothing
end

function calculate_overlap(A::HBM,m::Int,n::Int)
    nlocal, nb, nη, ns, nk,lg = A.nlocal, A.nb,A.nη, A.ns, A.latt.nk,A.lg
    A.Λ .= 0.0+0.0im
    tmpΛ = reshape(A.Λ,ns,nη,nb*nk,ns,nη,nb*nk)

    λkp = zeros(ComplexF64,nb*nk,nb*nk)
    ur = zeros(ComplexF64,nlocal*lg^2,nb*nk)
    for iη in 1:A.nη 
        ul = reshape(view(A.Uk,:,:,iη,:),:,nb*nk) 
        ur .= reshape( circshift(reshape(ul,nlocal,lg,lg,nb*nk),(0,-m,-n,0)) ,nlocal*lg^2,nb*nk)
        λkp .= ul' * ur
        for is in 1:A.ns 
            tmpΛ[is,iη,:,is,iη,:] = λkp 
        end 
    end
    return nothing
end