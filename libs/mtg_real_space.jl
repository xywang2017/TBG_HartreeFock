include("bmLL.jl")

mutable struct Coords 
    z::Vector{ComplexF64}
    x1::Vector{Float64}
    x2::Vector{Float64}
    lr::Int # number of real space points along one direction 
    nr::Int # total number of points  

    Coords() = new()
end 

function initCoords(A::Coords,params::Params,q::Int;lr::Int = 1)
    A.lr = lr 
    
    # A.x2 = collect(0:(q*lr-1)) ./ (lr)
    # A.x1 = collect(0:(lr-1)) ./ (lr)

    lmax = (lr-1) #/2
    A.x1 = collect(-lmax:lmax) ./ lr 
    lmax = (q*lr-1)/2
    A.x2 = collect(-lmax:lmax) ./ lr 

    A.z = (reshape(A.x1,:,1)*params.a1 .+ reshape(A.x2,1,:)*params.a2)[:]
    A.nr = length(A.z)
end

mutable struct MTG 
    # String for valley 
    _valley::String

    coord::Coords 

    params::Params 

    η::Int # valley degree of freedom  
    nl::Int # number of layers
    nq::Int # number of points (linear) in [0,1/q)

    p::Int 
    q::Int 
    qϕ::ComplexF64  # magnetic translation vector: ϕ/ϕ0 * 2π/L2 * ê_2
    nLL::Int   # 0,1,2... 
    nγ::Int    # 2, particle-hole symmetry 
    nH::Int    # nLL*nγ - 1

    lB::Float64  # magnetic length in absolute units 

    W::Array{ComplexF64,3}  # sublattice, nr, l γ n  η k

    fname::String # name to dump the MTG structure into 
    MTG() = new()
end

function constructMTG(bm::bmLL;lr::Int=1,fname::String="")
    A = MTG()
    A._valley = bm._valley
    A.η = isequal(A._valley,"K") ? 1 : -1
    A.p = bm.p 
    A.q = bm.q 
    A.qϕ = bm.qϕ 
    A.nLL = bm.nLL 
    A.nγ = bm.nγ 
    A.nH = bm.nH
    A.lB = bm.lB 
    A.params = bm.params 
    A.fname = fname
    A.nl = 2 # number of layers
    A.nq = bm.nq

    A.coord = Coords() 
    initCoords(A.coord,A.params,A.q; lr=lr)

    A.W = zeros(ComplexF64,2,A.nH*A.p*A.nl*A.q*A.nq^2,A.coord.nr)  # first 2 is sublattice A, B basis

    for iz in 1:size(A.W,3)
        constructMTGRealSpaceWavefunctions(iz,A,bm)
    end

    if !isempty(fname)
        jldopen(fname,"w") do file 
            file["MTG"] = A
        end
    end
    return A
end


function constructMTGRealSpaceWavefunctions(iz::Int,A::MTG,bm::bmLL)
    Wz = zeros(ComplexF64,2,A.nH,A.p,A.nl,A.q,A.nq,A.nq) # first 2 is sublattice basis
    _z = A.coord.z[iz]
    svec = collect(-20:20) # range of s to consider in generating MTG eigenstates from LL wavefunctions 
    for ik2 in 1:A.nq, l in 1:A.nl, r in 1:A.p, iH in 1:A.nH 
        Kr = (l==1) ? A.η*A.params.Kb : A.η*A.params.Kt
        n,γ = inγ(iH)
        p2 = bm.latt.k2[ik2+(r-1)*A.nq]
        for s in svec 
            _k = (bm.latt.k2[ik2+(r-1)*A.nq] - s*A.p/A.q ) * A.params.g2
            expfactor = exp(1im * 2π * s * (-p2*projector_para(A.params.a1,A.params.a2)/abs(A.params.a2)) ) * 
                            exp(1im *s^2/2 *projector_para(A.qϕ,A.params.a1)*abs(A.params.a1)) * 
                            exp(-1im * s * projector_norm(Kr,A.params.a2)*projector_norm(A.params.a1,A.params.a2)) 
            Φz = Ψ(_z,A.η,n,γ,l,_k,A.params,A.lB)
            for ik1 in 1:A.nq, iq in 1:A.q
                k1 = bm.latt.k2[ik1+(iq-1)*A.nq]
                Wz[:,iH,r,l,iq,ik1,ik2] += expfactor * exp(1im*2π*k1*s) * Φz 
            end
        end
    end
    A.W[:,:,iz] = reshape(Wz,2,:)
    return nothing 
end


function Ψ(z::ComplexF64,η::Int,n::Int,γ::Int,l::Int,k::ComplexF64,params::Params,lB::Float64)
    # η = ± 1 : valley index
    x,y = projector_norm(z,params.a2), projector_para(z,params.a2)
    Kl = (l==1) ? (3-2η)*params.Kb : (3-2η)*params.Kt
    xtilde = x/lB + projector_para(k-Kl,params.a2)*lB
    ψ = ComplexF64[0;0]
    θstrain = angle(params.a2)-π/2
    if η == 1 
        if n == 0
            ψ .= ComplexF64[0; Φ(n,xtilde)]
        else
            ψ .= ComplexF64[-1im*γ*exp(-1im*θstrain)*Φ(n-1,xtilde); Φ(n,xtilde)]/sqrt(2)
        end
    elseif η == -1 
        if n == 0
            ψ .= ComplexF64[Φ(n,xtilde); 0]
        else
            ψ .= ComplexF64[Φ(n,xtilde);1im*γ*exp(-1im*θstrain)*Φ(n-1,xtilde)]/sqrt(2)
        end
    else
        println("errow with η value")
    end

    ψ .= ψ * exp(1im*projector_norm(Kl,params.a2)*x) * exp(1im*projector_para(k,params.a2)*y)
    return ψ
end

function Φ(n::Int,x::Float64)
    if n < 150
        ψ=1/π^(1/4)*exp(-n/2*log(2)-x^2/2-0.5*sf_lnfact(n))*hermiteh(n,x)
    else
        # ψ=2^(-n/2)/(π^(1/4)*factorial(big(n))^(1/2))*exp(-x^2/2)*hermiteh(big(n),x)
        ψ = 0.0
    end
    return (abs(ψ)<1e-16) ? 0.0 : Float64(ψ)
end
