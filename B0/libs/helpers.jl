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