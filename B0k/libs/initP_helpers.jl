using Random 

function init_P(hf::HartreeFock; _Init::String="Random",
                    P0::Array{ComplexF64,3}=ones(ComplexF64,1,1,1))
    """
        Initialization of density matrix 
    """
    if isequal(_Init,"random") # random occupation of BM bands
        init_P_random(hf)
    elseif isequal(_Init,"flavor")
        init_P_flavor_polarization(hf)
    else
        hf.P .= P0 
    end
    init_P_random_rotation(hf;α=1.0)
    println("Initial filling is: ", real( 8*sum([tr(hf.P[:,:,ik]+0.5I) for ik in 1:size(hf.P,3)])/(size(hf.P,3)*size(hf.P,1))-4 ) )
    
    return nothing
end

function init_P_random_rotation(hf::HartreeFock;α::Float64=0.2)
    vecs = zeros(ComplexF64,hf.nt*hf.latt.nk,hf.nt*hf.latt.nk)
    vecs .= eigvecs(Hermitian(rand(ComplexF64,size(vecs))))
    hf.P .= (1-α)* hf.P .+ α * vecs' * hf.P * vecs
    return nothing
end

function init_P_random(hf::HartreeFock)
    # random initialization of density matrix for all flavors 
    νmax = round(Int,(hf.ν+4)/8 * size(hf.P,1))
    n_total = size(hf.P,1)
    states_to_populate = randperm(n_total)[1:νmax]
    for ip in states_to_populate
        hf.P[ip,ip] = 1.0 
    end
    hf.P .= hf.P - 0.5I
    return nothing
end

function init_P_flavor_polarization(hf::HartreeFock)
    return nothing
end