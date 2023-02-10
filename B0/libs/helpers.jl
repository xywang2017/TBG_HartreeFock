function check_Hermitian(H::Matrix{ComplexF64})
    err = norm(H - H')
    if err > 1e-6
        println("Error with Hermiticity: ",err)
    end
    return nothing
end

function check_Unitary(H::Matrix{ComplexF64})
    err = norm(H*H'-I)
    if err > 1e-6
        println("Error with Unitarity: ",err)
    end
    return nothing
end

function find_chemicalpotential(energies::Vector{Float64},ν::Float64)
    E = sort(energies)
    νs = collect(1:length(E)) ./ length(E)
    i = 1
    while ν > νs[i]
        i = i+1 
    end   
    # μ = (E[i-1]+E[i])/2
    μ = (E[i+1]+E[i])/2
    return μ
end