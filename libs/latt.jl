using LinearAlgebra

mutable struct Lattice 
    k1::Vector{Float64}
    k2::Vector{Float64}
    nk::Int # total k points in the mesh
    kvec::Vector{ComplexF64} # k1+ik2
    flag_inv::Bool # if true, the kvec is inverse symmetric

    Lattice() = new()
end

function constructLattice(latt::Lattice,params::Params;lk::Int=12)
    # even grid, does not go through Γ point
    # odd grid, does not go through Γ point
    lmax = (lk%2==0) ? (lk÷2-0.5) : ((lk-1)÷2)
    # latt.k1 = collect((-lmax):lmax) ./ lk 
    # latt.k2 = collect((-lmax):lmax) ./ lk 
    latt.k1 = collect(0:(lk-1)) ./ lk 
    latt.k2 = collect(0:(lk-1)) ./ lk 
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