include("BM_mod.jl")

mutable struct HartreeFock
    ns::Int 
    nb::Int 
    nη::Int 
    nt::Int

    fname::String # filename to load overlap matrix calculated from BM

    params::Params 
    latt::Lattice

    V0::Float64 # Coulomb unit
    ν::Float64 # filling fraction -4 to 4 
    P::Array{ComplexF64,3}  # one particle density matrix
    ϵk::Matrix{Float64} # eigenvalues of HF renormalized band dispersions nfl x ns x nk

    Λ::Array{ComplexF64,2}
    H0::Array{ComplexF64,3}
    H::Array{ComplexF64,3} # running Hamiltonian for any given k; nfl*ns x nfl*ns x nk
    precision::Float64  # iteration stopping point

    HartreeFock() = new()
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

@inline function V(q::ComplexF64,Lm::Float64) ::Float64
    res = 1e-6
    ϵr = 10.0
    return ( abs(q) < res ) ? 0 : 2π/(ϵr*abs(q))*tanh(abs(q)*Lm/2)
end

function run_HartreeFock(hf::HartreeFock,params::Params,latt::Lattice,fname::String;
            ν::Float64=0.0,savename::String="placeholder.txt",_Init::String="flavor")
    hf.ns, hf.nη, hf.nb, hf.nt = 2,2,2,8
    hf.params = params 
    hf.latt = latt

    hf.ν = ν 
    hf.precision = 1e-5

    hf.fname = fname

    hf.P = zeros(ComplexF64,hf.nt,hf.nt,latt.nk)
    hf.H = zeros(ComplexF64,size(hf.P))
    hf.ϵk = zeros(Float64,hf.nt,latt.nk)

    # Coulomb scale in meV units 
    hf.V0 = CoulombUnit(hf.params)

    init_P(hf,_Init=_Init)
    hf.Λ = zeros(ComplexF64,hf.nt*latt.nk,hf.nt*latt.nk)
    
    BM_info(hf)

    # Hartree-Fock iterations
    norm_convergence = 10.0 
    iter = 1
    iter_err = Float64[]
    iter_energy = Float64[]
    while norm_convergence > hf.precision
        println("Iter: ",iter)
        hf.H .= hf.H0 * 0.0
        add_Hartree(hf;β=1.0)
        add_Fock(hf;β=1.0)
        # @time add_Fock_vectorize(hf;β=1.0)
        Etot = compute_HF_energy(hf.H .- hf.H0,hf.H0,hf.P)

        #Δ is a projector to make it closed shell
        if norm_convergence <1e-4
            Δ = 0.0 
        else 
            Δ = 0.0
        end
        norm_convergence = update_P(hf;Δ=Δ,α=0.4)

        println("Running HF energy: ",Etot)
        println("Running norm convergence: ",norm_convergence)
        push!(iter_energy,Etot)
        push!(iter_err,norm_convergence)
        iter +=1
    end

    return nothing
end


function BM_info(hf::HartreeFock)
    hf.H0 = zeros(ComplexF64,size(hf.H))
    hbm = zeros(ComplexF64,hf.nη*hf.nb,hf.latt.nk)
    tmp_H = reshape(hf.H0,hf.ns,hf.nη*hf.nb,hf.ns,hf.nη*hf.nb,hf.latt.nk)
    
    jldopen(hf.fname,"r") do file
        hbm .= reshape(file["E"],hf.nη*hf.nb,hf.latt.nk)
        for ifl in 1:hf.nη*hf.nb, is in 1:hf.ns
            tmp_H[is,ifl,is,ifl,:] = hbm[ifl,:]
        end
    end
    return nothing
end

function add_Hartree(hf::HartreeFock;β::Float64=1.0)
    Gs = load(hf.fname,"Gs")
    lG = load(hf.fname,"lG")
    Glabels = (-(lG-1)÷2):((lG-1)÷2)
    Lm = sqrt(abs(hf.params.a1)*abs(hf.params.a2))

    tmp_Λ = reshape(hf.Λ,hf.nt,hf.latt.nk,hf.nt,hf.latt.nk)
    tmp_Λ1 = reshape(hf.Λ,hf.ns,hf.nη*hf.nb*hf.latt.nk,hf.ns,hf.nη*hf.nb*hf.latt.nk)
    
    for ig in 1:lG^2
        m,n = Glabels[(ig-1)%lG+1],Glabels[(ig-1)÷lG+1]
        jldopen(hf.fname,"r") do file 
            for is in 1:hf.ns 
                tmp_Λ1[is,:,is,:] = file["$(m)_$(n)"]
            end
        end
        trPG = 0.0+0.0im
        for ik in 1:hf.latt.nk 
            trPG += tr(view(hf.P,:,:,ik)*conj(view(tmp_Λ,:,ik,:,ik)))
        end
        for ik in 1:size(hf.H,3) 
            hf.H[:,:,ik] .+= ( β/hf.latt.nk*hf.V0*V(Gs[ig],Lm) * trPG) * view(tmp_Λ,:,ik,:,ik)
        end
    end
    return nothing
end

function add_Fock(hf::HartreeFock;β::Float64=1.0)
    Gs = load(hf.fname,"Gs")
    lG = load(hf.fname,"lG")
    Glabels = (-(lG-1)÷2):((lG-1)÷2)
    Lm = sqrt(abs(hf.params.a1)*abs(hf.params.a2))

    tmp_Λ = reshape(hf.Λ,hf.nt,hf.latt.nk,hf.nt,hf.latt.nk)
    tmp_Λ1 = reshape(hf.Λ,hf.ns,hf.nη*hf.nb*hf.latt.nk,hf.ns,hf.nη*hf.nb*hf.latt.nk)

    kvec = reshape(hf.latt.kvec,:) 
    tmp_Fock = zeros(ComplexF64,hf.nt,hf.nt)

    for ig in 1:lG^2 
        m,n = Glabels[(ig-1)%lG+1],Glabels[(ig-1)÷lG+1]
        jldopen(hf.fname,"r") do file 
            for is in 1:hf.ns 
                tmp_Λ1[is,:,is,:] = file["$(m)_$(n)"]
            end
        end
        for ik in 1:hf.latt.nk
            tmp_Fock .= 0.0 + 0.0im
            for ip in 1:hf.latt.nk
                tmp_Fock .+= ( β*hf.V0*V(kvec[ip]-kvec[ik]+Gs[ig],Lm) /hf.latt.nk) * 
                            ( view(tmp_Λ,:,ik,:,ip)*transpose(view(hf.P,:,:,ip))*view(tmp_Λ,:,ik,:,ip)' )
            end
            hf.H[:,:,ik] .-= tmp_Fock
        end
    end
    return nothing
end

function add_Fock_vectorize(hf::HartreeFock;β::Float64=1.0)
    Gs = load(hf.fname,"Gs")
    lG = load(hf.fname,"lG")
    Glabels = (-(lG-1)÷2):((lG-1)÷2)
    Lm = sqrt(abs(hf.params.a1)*abs(hf.params.a2))

    tmp_Λ = reshape(hf.Λ,hf.ns,hf.nη*hf.nb*hf.latt.nk,hf.ns,hf.nη*hf.nb*hf.latt.nk)
    kvec = reshape(hf.latt.kvec,:) 

    Λ_perm = zeros(ComplexF64,hf.nt*hf.latt.nk,hf.nt*hf.latt.nk)

    Vvals = zeros(Float64,hf.latt.nk,hf.latt.nk)
    λl = zeros(ComplexF64,1,hf.nt,hf.latt.nk,hf.nt,1,hf.latt.nk)
    λr = zeros(ComplexF64,hf.nt,1,hf.latt.nk,1,hf.nt,hf.latt.nk)
    vv = zeros(Float64,1,1,hf.latt.nk,1,1,hf.latt.nk)
    pp = zeros(ComplexF64,hf.nt,hf.nt,hf.latt.nk,1,1,1)
    
    for ig in 1:lG^2 
        m,n = Glabels[(ig-1)%lG+1],Glabels[(ig-1)÷lG+1]
        
        jldopen(hf.fname,"r") do file 
            for is in 1:hf.ns 
                tmp_Λ[is,:,is,:] = file["$(m)_$(n)"]
            end
        end
        Λ_perm .= transpose(hf.Λ)
        Vvals .= (β*hf.V0/hf.latt.nk) * V.(reshape(kvec,:,1) .- reshape(kvec,1,:) .+ Gs[ig],Lm)
        λl .= reshape(Λ_perm,1,hf.nt,hf.latt.nk,hf.nt,1,hf.latt.nk)
        λr .= conj.(reshape(Λ_perm,hf.nt,1,hf.latt.nk,1,hf.nt,hf.latt.nk))
        vv .= reshape(Vvals,1,1,hf.latt.nk,1,1,hf.latt.nk)
        pp .= reshape(hf.P,hf.nt,hf.nt,hf.latt.nk,1,1,1)

        for i3 in 1:hf.latt.nk, i2 in 1:hf.nt, i1 in 1:hf.nt 
            hf.H .-=  view(vv,1,1,i3,:,:,:) .*
                        view(λl,1,i2,i3,:,:,:).*
                        view(pp,i1,i2,i3,1,1,1).*
                        view(λr,i1,1,i3,:,:,:)
        end
    end
    return nothing
end

function update_P(hf::HartreeFock;Δ::Float64=0.0,α::Float64=0.3)
    """
        Diagonalize Hamiltonian for every k; use ν to keep the lowest N particle states;
        update P 
    """
    νnorm = round(Int,(hf.ν+4)/8 * size(hf.H,1)*size(hf.H,3))  # total number of occupied states 
    vecs = similar(hf.H)
    for ik in 1:size(hf.H,3)
        check_Hermitian(hf.H[:,:,ik])
        hf.ϵk[:,ik],vecs[:,:,ik] = eigen(Hermitian(view(hf.H,:,:,ik)-Δ*(conj.(view(hf.P,:,:,ik))+0.5*I)) )
    end

    iϵ_sorted = sortperm(hf.ϵk[:])
    iϵ_occupied = iϵ_sorted[1:νnorm]
    iband_occupied = (iϵ_occupied .-1) .% size(hf.ϵk,1) .+1
    ik_occupied = (iϵ_occupied .-1) .÷ size(hf.ϵk,1) .+1

    P_new = similar(hf.P)

    for ik in 1:size(hf.P,3)
        occupied_vecs = vecs[:,iband_occupied[ik_occupied.==ik],ik]
        P_new[:,:,ik] = conj(occupied_vecs)*transpose(occupied_vecs) - 0.5*I
    end

    norm_convergence = norm(P_new .- hf.P) ./ norm(P_new)
    hf.P .= α*P_new + (1-α)*hf.P
    
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
        println("Error with Hermitian Hamiltonian: ",err)
    end
    return nothing
end

include("initP_helpers.jl")