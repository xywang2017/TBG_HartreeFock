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
    σzτz::Matrix{Float64} # eigenvalues of σzτz in every HF state
    μ::Float64 # running Hartree-Fock chemical potential 
    Δ::Vector{Float64} # spin-valley-band mixing order parameter s_iη_jn_k
    Δstr::Vector{String} # string

    Λ::Array{ComplexF64,2}
    H0::Array{ComplexF64,3}
    Σz::Array{ComplexF64,3} # σzτz operator in the BM basis
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
    ϵr = 15.0
    return ( abs(q) < res ) ? 0 : 2π/(ϵr*abs(q))*tanh(abs(q)*4*Lm/2)
end

function run_HartreeFock(hf::HartreeFock,params::Params,latt::Lattice,fname::String;
            ν::Float64=0.0,savename::String="placeholder.txt",_Init::String="flavor")
    
    hf.params = params 
    hf.latt = latt

    hf.ν = ν 
    hf.precision = 1e-5

    hf.fname = fname
    jldopen(hf.fname,"r") do file 
        hf.ns, hf.nη, hf.nb = file["ns"],file["nη"],file["nb"]
        hf.nt = hf.ns*hf.nη*hf.nb
    end

    hf.P = zeros(ComplexF64,hf.nt,hf.nt,hf.latt.nk)
    hf.H = zeros(ComplexF64,size(hf.P))
    hf.Σz = zeros(ComplexF64,size(hf.P))
    hf.ϵk = zeros(Float64,hf.nt,latt.nk)
    hf.σzτz = zeros(Float64,hf.nt,latt.nk)

    # Coulomb scale in meV units 
    hf.V0 = CoulombUnit(hf.params)

    init_P(hf,_Init=_Init)
    hf.Λ = zeros(ComplexF64,hf.nt*latt.nk,hf.nt*latt.nk)
    
    BM_info(hf)

    # order parameters 
    γs = ["γ0","γ1","γ2","γ3"]
    σs = ["s0","s1","s2","s3"]
    τs = ["τ0","τ1","τ2","τ3"]
    hf.Δstr = [γs[i]*τs[j]*σs[k] for i in 1:4 for j in 1:4 for k in 1:4]
    hf.Δ = zeros(Float64,size(hf.Δstr))

    # Hartree-Fock iterations
    norm_convergence = 10.0 
    iter = 0
    iter_err = Float64[]
    iter_energy = Float64[]
    iter_oda = Float64[]
    while norm_convergence > hf.precision
        @time begin 
            
            hf.H .= hf.H0 * 1.0
            add_HartreeFock(hf;β=1.0)
            # add_Hartree(hf;β=1.0)
            # add_Fock(hf;β=1.0)
            # @time add_Fock_vectorize(hf;β=1.0)
            Etot = compute_HF_energy(hf.H .- hf.H0,hf.H0,hf.P)

            #Δ is a projector to make it closed shell -- incompatible with ODA
            Δ = 0.0
            norm_convergence,λ = update_P(hf;Δ=Δ)
            push!(iter_energy,Etot)
            push!(iter_err,norm_convergence)
            push!(iter_oda,λ)
            iter +=1
            if (mod(iter,10) == 1 )|| norm_convergence < hf.precision
                jldopen(savename,"w") do file 
                    file["hf"] = hf
                    file["iter_energy"] = iter_energy
                    file["iter_err"] = iter_err 
                    file["iter_oda"] = iter_oda
                end
            end

            if iter >300 || λ < 1e-3
                break 
            end

            println("Iter: ",iter)
            println("Running HF energy: ",Etot)
            println("Running norm convergence: ",norm_convergence)
            println("ODA parameter λ: ",λ)
        end
    end

    return nothing
end


function BM_info(hf::HartreeFock)
    hf.H0 = zeros(ComplexF64,size(hf.H))

    jldopen(hf.fname,"r") do file
        hf.H0 .= file["E"]
        hf.Σz .= file["Σz"]
    end
    return nothing
end

function add_Hartree(hf::HartreeFock;β::Float64=1.0)
    Gs = load(hf.fname,"Gs")
    lG = load(hf.fname,"lG")
    Glabels = (-(lG-1)÷2):((lG-1)÷2)
    Lm = sqrt(abs(hf.params.a1)*abs(hf.params.a2))

    tmp_Λ = reshape(hf.Λ,hf.nt,hf.latt.nk,hf.nt,hf.latt.nk)
    
    for ig in 1:lG^2
        m,n = Glabels[(ig-1)%lG+1],Glabels[(ig-1)÷lG+1]
        jldopen(hf.fname,"r") do file 
            hf.Λ .= file["$(m)_$(n)"]
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

    kvec = reshape(hf.latt.kvec,:) 
    tmp_Fock = zeros(ComplexF64,hf.nt,hf.nt)

    for ig in 1:lG^2 
        m,n = Glabels[(ig-1)%lG+1],Glabels[(ig-1)÷lG+1]
        jldopen(hf.fname,"r") do file 
            hf.Λ .= file["$(m)_$(n)"]
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

function add_HartreeFock(hf::HartreeFock;β::Float64=1.0)
    Gs = load(hf.fname,"Gs")
    lG = load(hf.fname,"lG")
    Glabels = (-(lG-1)÷2):((lG-1)÷2)
    Lm = sqrt(abs(hf.params.a1)*abs(hf.params.a2))

    tmp_Λ = reshape(hf.Λ,hf.nt,hf.latt.nk,hf.nt,hf.latt.nk)

    kvec = reshape(hf.latt.kvec,:) 
    tmp_Fock = zeros(ComplexF64,hf.nt,hf.nt)

    for ig in 1:lG^2 
        m,n = Glabels[(ig-1)%lG+1],Glabels[(ig-1)÷lG+1]
        G = Gs[ig]
        jldopen(hf.fname,"r") do file 
            hf.Λ .= file["$(m)_$(n)"]
        end

        G0 = abs(3*hf.params.g1+3*hf.params.g2)*1.00001
        if abs(G) <G0*cos(pi/6)/abs(cos(mod(angle(G),pi/3)-pi/6)) # this leads to a shell expansion up to 3g1+3g2
            ## Hartree
            trPG = 0.0+0.0im
            for ik in 1:hf.latt.nk 
                trPG += tr(view(hf.P,:,:,ik)*conj(view(tmp_Λ,:,ik,:,ik)))
            end
            for ik in 1:size(hf.H,3) 
                hf.H[:,:,ik] .+= ( β/hf.latt.nk*hf.V0*V(Gs[ig],Lm) * trPG) * view(tmp_Λ,:,ik,:,ik)
            end

            ## Fock
            for ik in 1:hf.latt.nk
                tmp_Fock .= 0.0 + 0.0im
                for ip in 1:hf.latt.nk
                    tmp_Fock .+= ( β*hf.V0*V(kvec[ip]-kvec[ik]+Gs[ig],Lm) /hf.latt.nk) * 
                                ( view(tmp_Λ,:,ik,:,ip)*transpose(view(hf.P,:,:,ip))*view(tmp_Λ,:,ik,:,ip)' )
                end
                hf.H[:,:,ik] .-= tmp_Fock
            end
        end
    end
    return nothing
end

function add_Fock_vectorize(hf::HartreeFock;β::Float64=1.0)
    Gs = load(hf.fname,"Gs")
    lG = load(hf.fname,"lG")
    Glabels = (-(lG-1)÷2):((lG-1)÷2)
    Lm = sqrt(abs(hf.params.a1)*abs(hf.params.a2))
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
            hf.Λ .= file["$(m)_$(n)"]
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

function update_P(hf::HartreeFock;Δ::Float64=0.0)
    """
        Diagonalize Hamiltonian for every k; use ν to keep the lowest N particle states;
        update P 
    """
    νnorm = round(Int,(hf.ν+4)/8 * size(hf.H,1)*size(hf.H,3))  # total number of occupied states 
    vecs = similar(hf.H)
    for ik in 1:size(hf.H,3)
        # check_Hermitian(hf.H[:,:,ik])
        # hf.ϵk[:,ik],vecs[:,:,ik] = eigen(Hermitian(view(hf.H,:,:,ik)-Δ*(conj.(view(hf.P,:,:,ik))+0.5*I)) )
        hf.ϵk[:,ik],vecs[:,:,ik] = eigen(Hermitian(view(hf.H,:,:,ik)))
        hf.σzτz[:,ik] = real( diag( view(vecs,:,:,ik)' * view(hf.Σz,:,:,ik) * view(vecs,:,:,ik) ) )
    end

    iϵ_sorted = sortperm(hf.ϵk[:])
    iϵ_occupied = iϵ_sorted[1:νnorm]
    iband_occupied = (iϵ_occupied .-1) .% size(hf.ϵk,1) .+1
    ik_occupied = (iϵ_occupied .-1) .÷ size(hf.ϵk,1) .+1

    hf.μ = find_chemicalpotential(hf.ϵk[:],(hf.ν+4)/8)
    hf.Δ .= calculate_valley_spin_band_order_parameters(hf)

    P_new = zeros(ComplexF64,size(hf.P))

    for ik in 1:size(hf.P,3)
        occupied_vecs = vecs[:,iband_occupied[ik_occupied.==ik],ik]
        P_new[:,:,ik] = conj(occupied_vecs)*transpose(occupied_vecs) - 0.5*I
    end

    norm_convergence = calculate_norm_convergence(P_new,hf.P)

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

function compute_HF_energy(H_HF::Array{ComplexF64,3},H0::Array{ComplexF64,3},P::Array{ComplexF64,3})
    Etot = 0.0 
    for ik in 1:size(H0,3) 
        Etot += tr(view(H_HF,:,:,ik)*transpose(view(P,:,:,ik)))/2 + tr(view(H0,:,:,ik)*transpose(view(P,:,:,ik)))
    end
    return real(Etot)/(size(H0,3))
end

function oda_parametrization(hf::HartreeFock,δP::Array{ComplexF64,3};β::Float64=1.0)
    # compute coefficients b λ + a λ^2/2

    Gs = load(hf.fname,"Gs")
    lG = load(hf.fname,"lG")
    Glabels = (-(lG-1)÷2):((lG-1)÷2)
    Lm = sqrt(abs(hf.params.a1)*abs(hf.params.a2))

    tmp_Λ = reshape(hf.Λ,hf.nt,hf.latt.nk,hf.nt,hf.latt.nk)
    kvec = reshape(hf.latt.kvec,:) 
    tmp_Fock = zeros(ComplexF64,hf.nt,hf.nt)

    # change of Hartree-Fock due to a small δP
    δH = zeros(ComplexF64,size(δP))
    for ig in 1:lG^2
        m,n = Glabels[(ig-1)%lG+1],Glabels[(ig-1)÷lG+1]
        G = Gs[ig]
        jldopen(hf.fname,"r") do file 
            hf.Λ .= file["$(m)_$(n)"]
        end
        G0 = abs(3*hf.params.g1+3*hf.params.g2)*1.00001
        if abs(G) <G0*cos(pi/6)/abs(cos(mod(angle(G),pi/3)-pi/6)) # this leads to a shell expansion up to 3g1+3g2
            trPG = 0.0+0.0im
            for ik in 1:hf.latt.nk 
                trPG += tr(view(δP,:,:,ik)*conj(view(tmp_Λ,:,ik,:,ik)))
            end
            for ik in 1:size(δH,3) 
                δH[:,:,ik] .+= ( β/hf.latt.nk*hf.V0*V(Gs[ig],Lm) * trPG) * view(tmp_Λ,:,ik,:,ik)
            end
            for ik in 1:hf.latt.nk
                tmp_Fock .= 0.0 + 0.0im
                for ip in 1:hf.latt.nk
                    tmp_Fock .+= ( β*hf.V0*V(kvec[ip]-kvec[ik]+Gs[ig],Lm) /hf.latt.nk) * 
                                ( view(tmp_Λ,:,ik,:,ip)*transpose(view(δP,:,:,ip))*view(tmp_Λ,:,ik,:,ip)' )
                end
                δH[:,:,ik] .-= tmp_Fock
            end
        end
    end

    # compute coefficients with δP 
    a, b = 0.0+0.0im , 0.0 +0.0im
    for ik in 1:hf.latt.nk 
        b += tr(transpose(view(δP,:,:,ik))*view(hf.H0,:,:,ik)) + 
             tr(transpose(view(δP,:,:,ik))*view(hf.H .- hf.H0,:,:,ik))/2 +
             tr(transpose(view(hf.P,:,:,ik))*view(δH,:,:,ik))/2
        a += tr(transpose(view(δP,:,:,ik))*view(δH,:,:,ik))
    end
    a = real(a)/size(δP,3)
    b = real(b)/size(δP,3)

    # if abs(imag(a))+abs(imag(b)) >1e-13 
    #     println(imag(a)," ",imag(b))
    # end
    λ, λ0 = 0.0 , -b/a
    # println("a= ",a," b= ",b," λ0=",λ0)
    if a>0 # convex and increasing with large λ 
        if λ0 <=0 
            λ = 0.0  # give it some kick..
        elseif λ0 <1 
            λ = λ0 
        else
            λ = 1.0 
        end
    else
        if λ0 <=0.5 
            λ = 1.0 
        else 
            λ = 0. # give it some kick..
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
    order_parameters = Float64[]
    Δ = zeros(Float64,size(hf.ϵk))
    for i in 1:4
        for j in 1:4 
            for k in 1:4
                for ik in 1:size(hf.ϵk,2)
                    F = eigen(Hermitian(view(hf.H,:,:,ik)))
                    Δ[:,ik] = real(diag(F.vectors'*kron(pauli_matrices[i],kron(pauli_matrices[j],pauli_matrices[k]))*F.vectors))
                end
                push!(order_parameters,sum(Δ[:][hf.ϵk[:].<= hf.μ])/length(hf.ϵk)*8)
            end
        end
    end
    return order_parameters
end


include("initP_helpers.jl")
include("helpers.jl")