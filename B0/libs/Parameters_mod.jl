using Parameters

@with_kw mutable struct Params
    vf::Float64 = 2680   # Liang number
    dθ::Float64 = 1.08π/180
    μ::Float64 = 0.0
    w0::Float64 = 77 ## AA tunneling
    w1::Float64 = 110 ## AB tunneling
    δ::Float64 = 0.0  # hBN alignment induced C2 symmetry breaking, only bottom layer

    # distance between Kt and Kb
    kb::Float64 = 8π/3*sin(dθ/2)

    # alternative choice of primary vectors - per Oskar & Jian choice
    g1::ComplexF64 = √3*kb
    g2::ComplexF64 = √3*kb*exp(1im * 2π/3)
    a1::ComplexF64 = 4π/(3kb)*exp(1im * π/6)
    a2::ComplexF64 = 4π/(3kb)*1im
    θ12::Float64 = π/3  # relative angle between a1 and a2 

    # coordinates for special points (Γ is at origin)
    Γ::ComplexF64 = 0.0im
    Kt::ComplexF64 = kb/2 * exp(1im * π/2)
    Kb::ComplexF64 = -Kt

    # if valley K prime, Kt and Kb are swapped

    # Tunneling matrix
    ω::ComplexF64 = exp(1im * 2π/3)   # - if valley K prime. + if valley K
    T0::Matrix{ComplexF64} = [[w0 w1];[w1 w0]]  # intra-unit cell, 
    T1::Matrix{ComplexF64} = [[w0 w1*conj(ω)];[w1*ω w0]]  # t -> b along +g2
    T2::Matrix{ComplexF64} = [[w0 w1*ω];[w1*conj(ω) w0]]  # t -> b along +g1

    # heterostrain
    ϵ::Float64 = 0.003
    φ::Float64 = 0.0*π/180
    ν::Float64 =  0.16
    ϵxx::Float64 = -ϵ * cos(φ)^2 + ν * ϵ * sin(φ)^2
    ϵyy::Float64 = ν * ϵ * cos(φ)^2 - ϵ * sin(φ)^2
    ϵxy::Float64 = (1+ν) * ϵ * cos(φ) * sin(φ)
    βg::Float64 = 3.14
    A::Vector{Float64} = (sqrt(3)*βg/2)*[ϵxx-ϵyy;-2ϵxy]
    Rφ::Matrix{Float64} = [cos(φ) -sin(φ);sin(φ) cos(φ)]
    S::Matrix{Float64} = Rφ' * [-ϵ 0; 0 ν*ϵ] * Rφ

    # alpha=0.5 means heterostrain, alpha=1/0 means strain on single layer
    α::Float64 = 0.5  
    # deformation potential strength 
    Da::Float64 = -4100.0 # meV
end

function initParamsWithStrain(params::Params)
    T = params.dθ/2 * Float64[0 -1; 1 0]
    G1, G2 = 4π/sqrt(3) * [0.0;-1], 4π/sqrt(3) * [sqrt(3)/2;0.5]
    tmp1 = (2T - params.S)*G1
    tmp2 = (2T - params.S)*G2
    params.g1 = tmp1[1] + 1im * tmp1[2]
    params.g2 = tmp2[1] + 1im * tmp2[2]
    
    area = abs(real(params.g1)*imag(params.g2)-imag(params.g1)*real(params.g2))
    params.a1 = 2π/area * (imag(params.g2)-1im*real(params.g2))
    params.a2 = 2π/area * (-imag(params.g1)+1im*real(params.g1))

    params.θ12 = angle(params.a2) - angle(params.a1)

    params.Kt = params.Kt + (params.A[1]+1im*params.A[2])*params.α - (params.S[1,1]+1im*params.S[2,1])*4π/3 * params.α
    params.Kb = params.Kb - (params.A[1]+1im*params.A[2])*(1-params.α) + (params.S[1,1]+1im*params.S[2,1])*4π/3 *(1-params.α)
    
     # this places Dirac cones at zone corners while still respecting particle hole symmetry
     params.Kt, params.Kb = params.Kt - params.g1/2, params.Kb + params.g1/2
    return nothing 

end