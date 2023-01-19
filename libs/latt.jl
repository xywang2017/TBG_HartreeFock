mutable struct Lattice 
    lk::Int # only lk=even is implemented
    listk::Matrix{Int} # indices of the k-grid in the first Moire Brillouin zone, i1,i2, (i2-1)lk+i1
    indΓ::Vector{Int}  # index of the Γ point in the mBZ, iΓ1,iΓ2,(iΓ2-1)lk+iΓ1
    kvec::Vector{ComplexF64} # values of kvectors in the first Moire Brilloun zone
    k1::Vector{Float64}
    k2::Vector{Float64}
    Lattice() = new()
end

function constructLattice(Latt::Lattice,params::Params;lk::Int=12)
    Latt.lk = lk
    itr = 1:Latt.lk
    Latt.listk = zeros(Int,3,lk^2)
    for i2 in itr, i1 in itr
        Latt.listk[3,(i2-1)*lk+i1] = (i2-1)*lk+i1
        Latt.listk[1,(i2-1)*lk+i1] = i1 
        Latt.listk[2,(i2-1)*lk+i1] = i2 
    end

    Latt.indΓ = Latt.listk[:,1]

    Latt.kvec = zeros(ComplexF64,lk^2)
    for ik in 1:lk^2
        Latt.kvec[ik] = (Latt.listk[1,ik] * params.g1 + Latt.listk[2,ik] * params.g2)/lk
    end
    # Latt.kvec .= Latt.kvec .- Latt.kvec[Latt.indΓ[3]] 
    
    Latt.k1 = collect((0.5-lk/2):(lk/2-0.5))/lk 
    Latt.k2 = collect((0.5-lk/2):(lk/2-0.5))/lk
    # Latt.k1 = collect(0:(lk-1))/lk 
    # Latt.k2 = collect(0:(lk-1))/lk
    Latt.kvec = (reshape(Latt.k1,:,1)*params.g1 .+ reshape(Latt.k2,1,:)*params.g2)[:]

    return nothing
end