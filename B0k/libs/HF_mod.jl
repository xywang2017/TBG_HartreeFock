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
    P::Array{ComplexF64,2}  # one particle density matrix
    ϵk::Vector{Float64} # eigenvalues of HF renormalized band dispersions nfl x ns x nk
    σzτz::Vector{Float64} # eigenvalues of σzτz in every HF state
    μ::Float64 # running Hartree-Fock chemical potential 
    Δ::Vector{Float64} # spin-valley-band mixing order parameter s_iη_jn_k
    Δstr::Vector{String} # string

    Λ::Array{ComplexF64,2}
    H0::Array{ComplexF64,2}
    Σz::Array{ComplexF64,2} # σzτz operator in the BM basis
    H::Array{ComplexF64,2} # running Hamiltonian for any given k; nfl*ns x nfl*ns x nk
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
    ϵr = 25.0
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

    hf.P = zeros(ComplexF64,hf.nt*hf.latt.nk,hf.nt*hf.latt.nk)
    hf.H = zeros(ComplexF64,size(hf.P))
    hf.Σz = zeros(ComplexF64,size(hf.P))
    hf.ϵk = zeros(Float64,hf.nt*latt.nk)
    hf.σzτz = zeros(Float64,hf.nt*latt.nk)

    # Coulomb scale in meV units 
    hf.V0 = CoulombUnit(hf.params)

    init_P(hf,_Init=_Init)
    hf.Λ = zeros(ComplexF64,hf.nt*latt.nk,hf.nt*latt.nk)
    
    BM_info(hf)

    # order parameters 
    ηs = ["η0","η1","η2","η3"]
    σs = ["s0","s1","s2","s3"]
    ns = ["n0","n1","n2","n3"]
    hf.Δstr = [ns[i]*ηs[j]*σs[k] for i in 1:4 for j in 1:4 for k in 1:4]
    hf.Δ = zeros(Float64,size(hf.Δstr))

    # Hartree-Fock iterations
    norm_convergence = 10.0 
    iter = 1
    iter_err = Float64[]
    iter_energy = Float64[]
    iter_oda = Float64[]
    while norm_convergence > hf.precision
        @time begin 
            println("Iter: ",iter)
            hf.H .= hf.H0 * 1.0
            add_HartreeFock(hf;β=1.0)
            Etot = compute_HF_energy(hf.H .- hf.H0,hf.H0,hf.P)

            #Δ is a projector to make it closed shell -- incompatible with ODA
            if norm_convergence <1e-4
                Δ = 0.0 
            else 
                Δ = 0.0
            end
            norm_convergence,λ = update_P(hf;Δ=Δ)
        end

        println("Running HF energy: ",Etot)
        println("Running norm convergence: ",norm_convergence)
        println("ODA parameter λ: ",λ)
        push!(iter_energy,Etot)
        push!(iter_err,norm_convergence)
        push!(iter_oda,λ)
        # if iter == 1 
        #     jldopen("typical_starting_point.jld2","w") do file 
        #         file["hf"] = hf 
        #     end
        # end
        iter +=1
        if (mod(iter,10) == 0 )|| norm_convergence < hf.precision
            jldopen(savename,"w") do file 
                file["hf"] = hf
                file["iter_energy"] = iter_energy
                file["iter_err"] = iter_err 
                file["iter_oda"] = iter_oda
            end
        end

        if iter >150 || λ < 1e-3
            break 
        end
    end

    return nothing
end


function BM_info(hf::HartreeFock)
    hf.H0 = zeros(ComplexF64,size(hf.H))
    tmpH0 = reshape(hf.H0,hf.nt,hf.latt.nk,hf.nt,hf.latt.nk)
    tmpΣz = reshape(hf.Σz,hf.nt,hf.latt.nk,hf.nt,hf.latt.nk)
    hbm = zeros(ComplexF64,hf.nt,hf.latt.nk)
    σz = zeros(ComplexF64,hf.nt,hf.nt,hf.latt.nk)

    jldopen(hf.fname,"r") do file
        hbm .= file["E"]
        σz .= file["Σz"]
        for ifl in 1:hf.nt, ik in 1:hf.latt.nk 
            tmpH0[ifl,ik,ifl,ik] = hbm[ifl,ik] 
            tmpΣz[ifl,ik,ifl,ik] = σz[ifl,ifl,ik]
        end
    end
    return nothing
end

function add_HartreeFock(hf::HartreeFock;β::Float64=1.0)
    Gs = load(hf.fname,"Gs")
    lG = load(hf.fname,"lG")
    Glabels = (-(lG-1)÷2):((lG-1)÷2)
    Lm = sqrt(abs(hf.params.a1)*abs(hf.params.a2))
    
    Λ1 = zeros(ComplexF64,hf.nt*hf.latt.nk,hf.nt*hf.latt.nk)
    Λ2 = zeros(ComplexF64,size(Λ1))
    tmpΛ1 = reshape(Λ1,hf.nt,hf.latt.lk,hf.latt.lk,hf.nt,hf.latt.lk,hf.latt.lk)
    tmpΛ2 = reshape(Λ2,size(tmpΛ1))
    δgs = [i+j*1im for i in -1:1 for j in -1:1]

    tmpH = reshape(hf.H,hf.nt,hf.latt.lk,hf.latt.lk,hf.nt,hf.latt.lk,hf.latt.lk)
    tmpP = reshape(hf.P,size(tmpH))

    fock = zeros(ComplexF64,hf.nt,hf.nt)
    λ1 = zeros(ComplexF64,hf.nt,hf.nt)
    λ2 = zeros(ComplexF64,size(λ1))
    p = zeros(ComplexF64,size(λ1))

    for ig in 1:lG^2 ,δg in δgs
        m,n = Glabels[(ig-1)%lG+1],Glabels[(ig-1)÷lG+1]
        m1,n1 = (m-Glabels[1]+real(δg))%lG + Glabels[1],(n-Glabels[1]+imag(δg))%lG + Glabels[1]
        jldopen(hf.fname,"r") do file 
            Λ1 .= file["$(m)_$(n)"]
            Λ2 .= file["$(m1)_$(n1)"]
        end
        for δk1 in (-(hf.latt.lk-1)):(hf.latt.lk-1),δk2 in (-(hf.latt.lk-1)):(hf.latt.lk-1)
            δp1,δp2 = δk1 - real(δg), δk2 - imag(δg)
            _δk = δk1*hf.params.g1 + δk2*hf.params.g2

            hartree = 0.0 + 0.0im 
            for ip2 in 1:hf.latt.lk, ip1 in 1:hf.latt.lk 
                ip1r, ip2r = ip1 -δp1, ip2 - δp2 
                if (ip1r in 1:hf.latt.lk && ip2r in 1:hf.latt.lk)
                    hartree += tr(conj(view(tmpΛ1,:,ip1r,ip2r,:,ip1,ip2))*view(tmpP,:,ip1,ip2,:,ip1r,ip2r))
                end
            end

            for ik2 in 1:hf.latt.lk, ik1 in 1:hf.latt.lk
                ik1r, ik2r = ik1 +δk1, ik2 + δk2 
                if (ik1r in 1:hf.latt.lk && ik2r in 1:hf.latt.lk)
                    vnum = (β*hf.V0/hf.latt.nk)*V(_δk+Gs[ig],Lm)
                    λ1 .= view(tmpΛ1,:,ik1,ik2,:,ik1r,ik2r)
                    fock .= 0.0 + 0.0im
                    for ip2 in 1:hf.latt.lk, ip1 in 1:hf.latt.lk 
                        ip1r, ip2r = ip1 -δp1, ip2 - δp2 
                        if (ip1r in 1:hf.latt.lk && ip2r in 1:hf.latt.lk)
                            λ2 .= view(tmpΛ2,:,ip1,ip2,:,ip1r,ip2r)
                            p .= view(tmpP,:,ip1,ip2,:,ik1r,ik2r)
                            fock .+= vnum * (λ1*transpose(p)*λ2')
                        end
                    end
                    tmpH[:,ik1,ik2,:,ik1r,ik2r] .+= hartree*λ1 - fock
                end
            end

        end 
    end
    return nothing
end

function update_P(hf::HartreeFock;Δ::Float64=0.0)
    """
        Diagonalize Hamiltonian for every k; use ν to keep the lowest N particle states;
        update P 
    """
    νnorm = round(Int,(hf.ν+4)/8 * size(hf.H,1))  # total number of occupied states 
    check_Hermitian(hf.H)
    hf.ϵk[:],vecs = eigen(Hermitian(hf.H))
    hf.σzτz .= real( diag( vecs' * hf.Σz * vecs ) )

    iϵ_sorted = sortperm(hf.ϵk)
    iϵ_occupied = iϵ_sorted[1:νnorm]

    hf.μ = find_chemicalpotential(hf.ϵk,(hf.ν+4)/8)
    hf.Δ .= calculate_valley_spin_band_order_parameters(hf)

    occupied_vecs = vecs[:,iϵ_occupied]
    P_new = conj(occupied_vecs)*transpose(occupied_vecs) - 0.5*I

    norm_convergence = calculate_norm_convergence(P_new,hf.P)

    λ = oda_parametrization(hf,P_new .- hf.P;β=1.0)
    # λ = 1.0 # often times oda_parameterization returns λ = 1.0, therefore not necessary
    norm_convergence = calculate_norm_convergence(λ*P_new + (1-λ)*hf.P,hf.P)
    hf.P .= λ*P_new + (1-λ)*hf.P
    return norm_convergence,λ
end

function calculate_norm_convergence(P2::Array{ComplexF64,2},P1::Array{ComplexF64,2})
    return norm(P1 .- P2) ./ norm(P2)
end

function compute_HF_energy(H_HF::Array{ComplexF64,2},H0::Array{ComplexF64,2},P::Array{ComplexF64,2})
    Etot = tr(H_HF*transpose(P))/2 + tr(H0*transpose(P))
    return real(Etot)/(size(H0,2))*8
end

function oda_parametrization(hf::HartreeFock,δP::Array{ComplexF64,2};β::Float64=1.0)
    # compute coefficients b λ + a λ^2/2
    Gs = load(hf.fname,"Gs")
    lG = load(hf.fname,"lG")
    Glabels = (-(lG-1)÷2):((lG-1)÷2)
    Lm = sqrt(abs(hf.params.a1)*abs(hf.params.a2))
    
    Λ1 = zeros(ComplexF64,hf.nt*hf.latt.nk,hf.nt*hf.latt.nk)
    Λ2 = zeros(ComplexF64,size(Λ1))
    tmpΛ1 = reshape(Λ1,hf.nt,hf.latt.lk,hf.latt.lk,hf.nt,hf.latt.lk,hf.latt.lk)
    tmpΛ2 = reshape(Λ2,size(tmpΛ1))
    δgs = [i+j*1im for i in -1:1 for j in -1:1]
    
    tmpP = reshape(δP,size(tmpH))
    # change of Hartree-Fock due to a small δP
    δH = zeros(ComplexF64,size(δP))
    tmpH = reshape(δH,hf.nt,hf.latt.lk,hf.latt.lk,hf.nt,hf.latt.lk,hf.latt.lk)

    fock = zeros(ComplexF64,hf.nt,hf.nt)
    λ1 = zeros(ComplexF64,hf.nt,hf.nt)
    λ2 = zeros(ComplexF64,size(λ1))
    p = zeros(ComplexF64,size(λ1))
    
    for ig in 1:lG^2 ,δg in δgs
        m,n = Glabels[(ig-1)%lG+1],Glabels[(ig-1)÷lG+1]
        m1,n1 = (m-Glabels[1]+real(δg))%lG + Glabels[1],(n-Glabels[1]+imag(δg))%lG + Glabels[1]
        jldopen(hf.fname,"r") do file 
            Λ1 .= file["$(m)_$(n)"]
            Λ2 .= file["$(m1)_$(n1)"]
        end
        for δk1 in (-(hf.latt.lk-1)):(hf.latt.lk-1),δk2 in (-(hf.latt.lk-1)):(hf.latt.lk-1)
            δp1,δp2 = δk1 - real(δg), δk2 - imag(δg)
            _δk = δk1*hf.params.g1 + δk2*hf.params.g2

            hartree = 0.0 + 0.0im 
            for ip2 in 1:hf.latt.lk, ip1 in 1:hf.latt.lk 
                ip1r, ip2r = ip1 -δp1, ip2 - δp2 
                if (ip1r in 1:hf.latt.lk && ip2r in 1:hf.latt.lk)
                    hartree += tr(conj(view(tmpΛ1,:,ip1r,ip2r,:,ip1,ip2))*view(tmpP,:,ip1,ip2,:,ip1r,ip2r))
                end
            end

            for ik2 in 1:hf.latt.lk, ik1 in 1:hf.latt.lk
                ik1r, ik2r = ik1 +δk1, ik2 + δk2 
                if (ik1r in 1:hf.latt.lk && ik2r in 1:hf.latt.lk)
                    vnum = (β*hf.V0/hf.latt.nk)*V(_δk+Gs[ig],Lm)
                    λ1 .= view(tmpΛ1,:,ik1,ik2,:,ik1r,ik2r)
                    fock .= 0.0 + 0.0im
                    for ip2 in 1:hf.latt.lk, ip1 in 1:hf.latt.lk 
                        ip1r, ip2r = ip1 -δp1, ip2 - δp2 
                        if (ip1r in 1:hf.latt.lk && ip2r in 1:hf.latt.lk)
                            λ2 .= view(tmpΛ2,:,ip1,ip2,:,ip1r,ip2r)
                            p .= view(tmpP,:,ip1,ip2,:,ik1r,ik2r)
                            fock .+= vnum * (λ1*transpose(p)*λ2')
                        end
                    end
                    tmpH[:,ik1,ik2,:,ik1r,ik2r] .+= hartree*λ1 - fock
                end
            end

        end 
    end

    # compute coefficients with δP 
    b = tr(transpose(δP)*hf.H0) + 
             tr(transpose(δP)*(hf.H .- hf.H0))/2 +
             tr(transpose(hf.P)*δH)/2
    a = tr(transpose(δP)*δH)
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
    F = eigen(Hermitian(hf.H))
    Ki = Array{ComplexF64,2}(I,hf.latt.nk,hf.latt.nk)
    for i in 1:4
        for j in 1:4 
            for k in 1:4
                Δ = real(diag(F.vectors'*kron(Ki,kron(pauli_matrices[i],kron(pauli_matrices[j],pauli_matrices[k])))*F.vectors))
                push!(order_parameters,sum(Δ[:][hf.ϵk.<= hf.μ])/length(hf.ϵk)*8)
            end
        end
    end
    return order_parameters
end


include("initP_helpers.jl")
include("helpers.jl")