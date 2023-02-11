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
    # valley x spin U(4) rotation --- otherwise above initializations do not access valley spin coherent states
    # init_P_valley_spin_rotation(hf;α=0.2)
    # init_P_valley_rotation(hf;α=0.2)
    init_P_random_rotation(hf;α=0.9)
    println("Initial filling is: ", real( 8*sum([tr(hf.P[:,:,ik]+0.5I) for ik in 1:size(hf.P,3)])/(size(hf.P,3)*size(hf.P,1))-4 ) )
    
    return nothing
end

function init_P_valley_spin_rotation(hf::HartreeFock;α::Float64=0.2)
    P0 = reshape(hf.P,hf.nη*hf.ns,hf.nb,hf.nη*hf.ns,hf.nb,hf.latt.nk)
    vecs = zeros(ComplexF64,hf.nη*hf.ns,hf.nη*hf.ns)
    for ik in 1:size(hf.P,3), iband in 1:size(P0,2)
        vecs .= eigvecs(Hermitian(rand(ComplexF64,hf.nη*hf.ns,hf.nη*hf.ns)))
        P0[:,iband,:,iband,ik] .= (1-α) *view(P0,:,iband,:,iband,ik) .+ α* vecs' * view(P0,:,iband,:,iband,ik) * vecs
    end
    return nothing
end

function init_P_valley_rotation(hf::HartreeFock;α::Float64=0.2)
    P0 = reshape(hf.P,hf.ns,hf.nb*hf.nη,hf.ns,hf.nb*hf.nη,hf.latt.nk)
    vecs = zeros(ComplexF64,hf.nb*hf.nη,hf.nb*hf.nη)
    for ik in 1:size(hf.P,3), is in 1:hf.ns
        vecs .= eigvecs(Hermitian(rand(ComplexF64,size(vecs))))
        P0[is,:,is,:,ik] .= (1-α) *view(P0,is,:,is,:,ik) .+ α* vecs' * view(P0,is,:,is,:,ik) * vecs
    end
    return nothing
end

function init_P_random_rotation(hf::HartreeFock;α::Float64=0.2)
    vecs = zeros(ComplexF64,hf.nt,hf.nt)
    for ik in 1:size(hf.P,3)
        vecs .= eigvecs(Hermitian(rand(ComplexF64,size(vecs))))
        hf.P[:,:,ik] .= (1-α)* view(hf.P,:,:,ik) .+ α * vecs' * view(hf.P,:,:,ik) * vecs
    end
    return nothing
end

function init_P_random(hf::HartreeFock)
    # random initialization of density matrix for all flavors 
    νmax = round(Int,(hf.ν+4)/8 * size(hf.P,1) * size(hf.P,3))
    n_total = size(hf.P,1) * size(hf.P,3)
    states_to_populate = randperm(n_total)[1:νmax]
    for ip in states_to_populate
        iq = (ip-1)%size(hf.P,1) +1 
        ik = (ip-1)÷size(hf.P,1) +1 
        hf.P[iq,iq,ik] = 1.0 
    end
    for ik in 1:size(hf.P,3)
        hf.P[:,:,ik] .= hf.P[:,:,ik] - 0.5I
    end
    return nothing
end

function init_P_flavor_polarization(hf::HartreeFock)
    # maximal population of the Chern states, and random sprinkling of the remaining states 
    νmax = round(Int,(hf.ν+4)/8 * size(hf.H,1) * size(hf.H,3))
    n_per_valley_spin = size(hf.H,1)*size(hf.H,3)÷(hf.nη*hf.ns)

    tmpP = reshape(hf.P,hf.nη*hf.ns,hf.nb,hf.nη*hf.ns,hf.nb,:)

    num_full_flavors_to_populate = νmax ÷ n_per_valley_spin 
    num_part_flavors_to_populate = (νmax % n_per_valley_spin ==0) ? 0 : 1
    
    # flavors_to_populate = randperm(hf.nη*hf.ns)[1:(num_full_flavors_to_populate+num_part_flavors_to_populate)]
    flavors_to_populate = collect(1:hf.nη*hf.ns)[1:(num_full_flavors_to_populate+num_part_flavors_to_populate)]
    # flavors_to_populate = [1;3]
    # println(num_part_flavors_to_populate)
    for iηs in flavors_to_populate[1:(end-num_part_flavors_to_populate)], iq in 1:size(tmpP,2)
        tmpP[iηs,iq,iηs,iq,:] .= 1.0 
    end
    if num_part_flavors_to_populate!=0
        iηs = flavors_to_populate[end]
        n_remain = νmax - (length(flavors_to_populate)-1)*n_per_valley_spin 
        states_to_populate = randperm(n_per_valley_spin)[1:n_remain]
        for ip in states_to_populate 
            iband = (ip-1)%(hf.nb) + 1 
            ik = (ip-1)÷(hf.nb) + 1
            tmpP[iηs,iband,iηs,iband,ik] = 1.0 
        end
    end
    
    for ik in 1:size(hf.P,3)
        hf.P[:,:,ik] .= hf.P[:,:,ik] - 0.5I
    end
    return nothing
end