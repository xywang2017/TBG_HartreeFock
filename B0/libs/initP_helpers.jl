using Random 

function init_P(hf::HartreeFock; _Init::String="Random",
                    P0::Array{ComplexF64,3}=ones(ComplexF64,1,1,1))
    """
        Initialization of density matrix 
    """
    if isequal(_Init,"random") # random occupation of BM bands
        init_P_random(hf)
    elseif isequal(_Init,"tivc")
        init_P_flavor_polarization(hf;flag="vp")
        init_P_ivc_rotation(hf)
    elseif isequal(_Init,"kivc")
        init_P_flavor_polarization(hf;flag="sublattice")
        init_P_ivc_rotation(hf)
    elseif isequal(_Init,"vp")
        init_P_flavor_polarization(hf;flag="vp")
    elseif isequal(_Init,"sp")
        init_P_flavor_polarization(hf;flag="sp")
    elseif isequal(_Init,"chern")
        init_P_flavor_polarization(hf;flag="chern")
    elseif isequal(_Init,"flavor")
        init_P_flavor_polarization(hf;flag="random")
    else
        hf.P .= P0 
    end
    # valley x spin U(4) rotation --- otherwise above initializations do not access valley spin coherent states
    # init_P_valley_spin_rotation(hf;α=0.2)
    # init_P_valley_rotation(hf;α=1.0)
    if isequal(_Init,"random")
        init_P_random_rotation(hf;α=1.0)
    end
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

function init_P_ivc_rotation(hf::HartreeFock)
    # rotate to tivc state from valley polarized state or sublattice polarized state(τx)
    P0 = reshape(hf.P,hf.ns,hf.nη,hf.nb,hf.ns,hf.nη,hf.nb,hf.latt.nk)
    mat = ComplexF64[1 1;1 -1]/sqrt(2)
    for ik in 1:size(hf.P,3), is in 1:hf.ns, ib in 1:hf.nb
        P0[is,:,ib,is,:,ib,ik] .= mat*view(P0,is,:,ib,is,:,ib,ik)*mat
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

function init_P_flavor_polarization(hf::HartreeFock;flag::String="vp")
    # maximal population of the Chern states, and random sprinkling of the remaining states 
    νmax = round(Int,(hf.ν+4)/8 * size(hf.H,1) * size(hf.H,3))
    n_per_flavor = size(hf.H,1)*size(hf.H,3)÷(hf.nt)

    num_full_flavors_to_populate = νmax ÷ n_per_flavor 
    num_part_flavors_to_populate = (νmax % n_per_flavor ==0) ? 0 : 1
    
    # 
    flavors = collect(reshape(1:hf.nt,hf.ns,hf.nη,hf.nb))
    if isequal(flag,"random")
        flavors_to_populate = randperm(hf.nt)[1:(num_full_flavors_to_populate+num_part_flavors_to_populate)]
    elseif isequal(flag,"vp")
        flavors_to_populate = permutedims(flavors,(1,3,2))[1:(num_full_flavors_to_populate+num_part_flavors_to_populate)]
    elseif isequal(flag,"sp")
        flavors_to_populate = permutedims(flavors,(3,2,1))[1:(num_full_flavors_to_populate+num_part_flavors_to_populate)]
    elseif isequal(flag,"chern")
        flavors_to_populate = permutedims(flavors,(1,2,3))[1:(num_full_flavors_to_populate+num_part_flavors_to_populate)]
    elseif isequal(flag,"sublattice")
        tmp = flavors[:,1,1]
        flavors[:,1,1] = flavors[:,1,2]
        flavors[:,1,2] = tmp
        flavors_to_populate = permutedims(flavors,(1,2,3))[1:(num_full_flavors_to_populate+num_part_flavors_to_populate)]
    end

    for ifl in flavors_to_populate[1:(end-num_part_flavors_to_populate)]
        hf.P[ifl,ifl,:] .= 1.0 
    end
    if num_part_flavors_to_populate!=0
        ifl = flavors_to_populate[end]
        n_remain = νmax - (length(flavors_to_populate)-1)*n_per_flavor 
        states_to_populate = randperm(n_per_flavor)[1:n_remain]
        for ip in states_to_populate 
            hf.P[ifl,ifl,ip] = 1.0 
        end
    end
    
    for ik in 1:size(hf.P,3)
        hf.P[:,:,ik] .= hf.P[:,:,ik] - 0.5I
    end
    return nothing
end