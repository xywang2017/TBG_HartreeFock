using LinearAlgebra
using Random
mutable struct Lattice 
    k1::Vector{Float64}
    k2::Vector{Float64}
    nk::Int # total k points in the mesh
    kvec::Vector{ComplexF64} # k1+ik2
    flag_inv::Bool # if true, the kvec is inverse symmetric

    Lattice() = new()
end

function constructLattice(latt::Lattice,params::Params;lk::Int=12)
    # lmax = (lk-1)/2
    # latt.k1 = collect((-lmax):lmax) ./ lk 
    # latt.k2 = collect((-lmax):lmax) ./ lk 
    latt.k1 = collect(0:(lk-1)) ./ lk
    latt.k2 = collect(0:(lk-1)) ./ lk 
    latt.kvec = (reshape(latt.k1,:,1)*params.g1 .+ reshape(latt.k2,1,:)*params.g2)[:]
    latt.nk = length(latt.kvec)
    latt.flag_inv = true
    return nothing
end

function constructLatticeIKS(latt::Lattice,params::Params;lk::Int=12,_valley::String="K",q0::Complex{Int}=0+0im)
    # even grid, does not go through Γ point
    # odd grid, does not go through Γ point
    # lmax = (lk-1)/2
    # latt.k1 = collect((-lmax):lmax) ./ lk 
    # latt.k2 = collect((-lmax):lmax) ./ lk 
    δk = rand(Float64) / lk 
    latt.k1 = collect(0:(lk-1)) ./ lk .+δk
    latt.k2 = collect(0:(lk-1)) ./ lk .+δk

    if isequal(_valley,"Kprime")
        latt.k1 .+= real(q0)/lk 
        latt.k2 .+= imag(q0)/lk
    end
    latt.kvec = (reshape(latt.k1,:,1)*params.g1 .+ reshape(latt.k2,1,:)*params.g2)[:]
    latt.nk = length(latt.kvec)
    latt.flag_inv = true
    return nothing
end


function initLatticeWithKvec(kvec::Vector{ComplexF64})
    latt = Lattice()
    latt.kvec = deepcopy(kvec)
    latt.k1 = real(latt.kvec)
    latt.k2 = imag(latt.kvec)
    latt.nk = length(latt.kvec)
    latt.flag_inv = false  # generally expect non-inversion-symmetric construct
    return latt
end