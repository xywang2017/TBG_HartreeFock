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
    if i < length(E)
        μ = (E[i+1]+E[i])/2
    else
        μ = E[i]
    end
    return μ
end



function find_lowest_energy_datafile(dir::String;test_str::String="null",_printinfo::Bool=false)
    metadata = ""
    if ! ispath(dir)
        # println("Hartree-Fock not run for the directory: ",dir)
    else 
        metadatas = String[]
        files = readdir(dir)
        for f in files 
            if occursin(test_str,f)
                push!(metadatas,joinpath(dir,f))
            end
        end
        if isempty(metadatas)
            # println("Hartree-Fock not run for the parameterization")
        else
            E = load(metadatas[1],"iter_energy")[end]
            metadata = metadatas[1]
            if length(metadatas)>1
                for i in 2:length(metadatas)
                    E0 = load(metadatas[i],"iter_energy")[end]
                    if E0<E 
                        E, metadata = E0, metadatas[i]
                    end
                end
            end
        end
    end
    if !isempty(metadata) && _printinfo
        println("HF energy: ",load(metadata,"iter_energy")[end])
        println("Convergence: ",load(metadata,"iter_err")[end])
    end
    return metadata
end
