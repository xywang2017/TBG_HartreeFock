
function init_P(hf::HartreeFock; _Init::String="BM",
    P0::Array{ComplexF64,3}=ones(ComplexF64,1,1,1),H0::Array{ComplexF64,3}=ones(ComplexF64,1,1,1))
    """
        Initialization of density matrix 
    """
    if isequal(_Init,"Random") # random occupation of BM bands
        init_P_random(hf)
    elseif isequal(_Init,"U Spectra") # filling of excitation spectra of strong coupling Hamiltonian at CNP
        init_P_strong_coupling(hf,H0=H0,P0=P0)
    elseif isequal(_Init,"Flavor U(4)")
        init_P_flavor_polarization(hf)
    elseif isequal(_Init,"Chern") # filling of Chern spectra of strong coupling Hamiltonian at CNP
        init_P_chern(hf,H0=H0)
    elseif isequal(_Init,"Sublattice")
        # init_P_sublattice(hf)
        init_P_sublattice_no_momentum(hf)
    else
        hf.P .= P0 
    end
    # valley x spin U(4) rotation --- otherwise above initializations do not access valley spin coherent states
    # init_P_valley_spin_roation(hf;α=0.2)
    init_P_random_rotation(hf;α=0.1)
    # init_P_valley_rotation(hf;α=0.2)
    println("Initial filling is: ", real( 8*sum([tr(hf.P[:,:,ik]+0.5I) for ik in 1:size(hf.P,3)])/(size(hf.P,3)*size(hf.P,1))-4 ) )
    
    return nothing
end

function init_P_valley_spin_rotation(hf::HartreeFock;α::Float64=0.2)
    P0 = reshape(hf.P,hf.nb*hf.q,hf.nη*hf.ns,hf.nb*hf.q,hf.nη*hf.ns,hf.q*hf.nq^2)
    vecs = zeros(ComplexF64,hf.nη*hf.ns,hf.nη*hf.ns)
    for ik in 1:size(hf.P,3), iband in 1:size(P0,1)
        vecs .= eigvecs(Hermitian(rand(ComplexF64,hf.nη*hf.ns,hf.nη*hf.ns)))
        P0[iband,:,iband,:,ik] .= (1-α)*view(P0,iband,:,iband,:,ik) + α* vecs' * view(P0,iband,:,iband,:,ik) * vecs
    end
    return nothing
end

function init_P_valley_rotation(hf::HartreeFock;α::Float64=0.2)
    P0 = reshape(hf.P,hf.nb*hf.q*hf.nη,hf.ns,hf.nb*hf.q*hf.nη,hf.ns,hf.q*hf.nq^2)
    vecs = zeros(ComplexF64,hf.nb*hf.q*hf.nη,hf.nb*hf.q*hf.nη)
    for ik in 1:size(hf.P,3), is in 1:hf.ns
        vecs .= eigvecs(Hermitian(rand(ComplexF64,size(vecs))))
        P0[:,is,:,is,ik] .= (1-α) *view(P0,:,is,:,is,ik) .+ α* vecs' * view(P0,:,is,:,is,ik) * vecs
    end
    return nothing
end

function init_P_random_rotation(hf::HartreeFock;α::Float64=0.2)
    vecs = zeros(ComplexF64,hf.nb*hf.q*hf.nη*hf.ns,hf.nb*hf.q*hf.nη*hf.ns)
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

function init_P_chern(hf::HartreeFock;H0::Array{ComplexF64,3}=ones(ComplexF64,1,1,1))
    # this function initializes the density matrix into a Chern state of the strong coupling spectrum 
    # first need to recreate the density matrix based on CNP valley-polarized Hamiltonian
    vecs = zeros(ComplexF64,size(H0,1),size(H0,2),size(H0,3)) 
    σzτz = zeros(Float64,size(H0,1),size(H0,3))
    idx_chern_minus = Complex{Int}[]
    idx_chern_plus = Complex{Int}[]
    for ik in 1:size(H0,3)
         vecs[:,:,ik] = eigvecs(Hermitian(view(H0,:,:,ik)))
         for iq in 1:size(H0,1)
            σzτz[iq,ik] = real(view(vecs,:,iq,ik)'*view(hf.Σz0,:,:,ik)*view(vecs,:,iq,ik))
         end
         idx = sortperm(σzτz[:,ik])
         for iq in 1:(4hf.q+4hf.p)
            push!(idx_chern_minus,idx[iq]+1im*ik)
         end
         for iq in (4hf.q+4hf.p+1):length(idx)
            push!(idx_chern_plus,idx[iq]+1im*ik)
         end
    end

    # maximal population of the Chern states, and random sprinkling of the remaining states 
    νmax = round(Int,(hf.ν+4)/8 * size(H0,1) * size(H0,3))
    states_to_populate = zeros(Complex{Int},νmax)
    if νmax < length(idx_chern_minus)
        states_to_populate .= shuffle(idx_chern_minus)[1:νmax]
    else
        states_to_populate[1:length(idx_chern_minus)] .= idx_chern_minus
        states_to_populate[(length(idx_chern_minus)+1):end] .= shuffle(idx_chern_plus)[1:(νmax-length(idx_chern_minus))]  
    end
    for ip in states_to_populate
        iband, ik = real(ip), imag(ip)
        hf.P[:,:,ik] .+= conj(view(vecs,:,iband,ik))*transpose(view(vecs,:,iband,ik))
    end

    for ik in 1:size(H0,3)
        hf.P[:,:,ik] .= hf.P[:,:,ik] - 0.5I
    end
    return nothing
end


function init_P_sublattice(hf::HartreeFock)
    # this function initializes the density matrix into a Chern state of the strong coupling spectrum 
    # first need to recreate the density matrix based on CNP valley-polarized Hamiltonian
    H0 = hf.Σz0
    vecs = zeros(ComplexF64,size(H0,1),size(H0,2),size(H0,3)) 
    σzτz = zeros(Float64,size(H0,1),size(H0,3))
    idx_chern_minus = Complex{Int}[]
    idx_chern_plus = Complex{Int}[]
    for ik in 1:size(H0,3)
         vecs[:,:,ik] = eigvecs(Hermitian(view(H0,:,:,ik)))
         for iq in 1:size(H0,1)
            σzτz[iq,ik] = real(view(vecs,:,iq,ik)'*view(hf.Σz0,:,:,ik)*view(vecs,:,iq,ik))
         end
         idx = sortperm(σzτz[:,ik])
         for iq in 1:(4hf.q+4hf.p)
            push!(idx_chern_minus,idx[iq]+1im*ik)
         end
         for iq in (4hf.q+4hf.p+1):length(idx)
            push!(idx_chern_plus,idx[iq]+1im*ik)
         end
    end

    # maximal population of the Chern states, and random sprinkling of the remaining states 
    νmax = round(Int,(hf.ν+4)/8 * size(H0,1) * size(H0,3))
    states_to_populate = zeros(Complex{Int},νmax)
    if νmax < length(idx_chern_minus)
        states_to_populate .= shuffle(idx_chern_minus)[1:νmax]
    else
        states_to_populate[1:length(idx_chern_minus)] .= idx_chern_minus
        states_to_populate[(length(idx_chern_minus)+1):end] .= shuffle(idx_chern_plus)[1:(νmax-length(idx_chern_minus))]  
    end
    for ip in states_to_populate
        iband, ik = real(ip), imag(ip)
        hf.P[:,:,ik] .+= conj(view(vecs,:,iband,ik))*transpose(view(vecs,:,iband,ik))
    end

    for ik in 1:size(hf.P,3)
        hf.P[:,:,ik] .= hf.P[:,:,ik] - 0.5I
    end
    return nothing
end


function init_P_sublattice_no_momentum(hf::HartreeFock)
    # this function initializes the density matrix into a Chern state of the strong coupling spectrum 
    # unlike init_P_sublattice, this function equally populates all k states, so only works for certain filling sequences
    H0 = hf.Σz0
    vecs = zeros(ComplexF64,size(H0,1),size(H0,2),size(H0,3)) 
    tmpH0 = reshape(H0,hf.nb*hf.q,hf.nη,hf.ns,hf.nb*hf.q,hf.nη,hf.ns,size(H0,3))
    tmpvecs = reshape(vecs,hf.nb*hf.q,hf.nη,hf.ns,hf.nb*hf.q,hf.nη,hf.ns,size(H0,3))
    σzτz = zeros(Float64,size(H0,1),size(H0,3))
    idx_chern_minus = Int64[]
    idx_chern_plus = Int64[]
    for ik in 1:size(H0,3)
        for iη in 1:hf.nη, is in 1:hf.ns
            tmpvecs[:,iη,is,:,iη,is,ik] = eigvecs(Hermitian(view(tmpH0,:,iη,is,:,iη,is,ik)))
        end
        for iq in 1:size(H0,1)
            σzτz[iq,ik] = real(view(vecs,:,iq,ik)'*view(H0,:,:,ik)*view(vecs,:,iq,ik))
        end
    end

    idx = sortperm(σzτz[:,1])
    for iq in 1:(4hf.q+4hf.p)
        push!(idx_chern_minus,idx[iq])
    end
    for iq in (4hf.q+4hf.p+1):length(idx)
        push!(idx_chern_plus,idx[iq])
    end

    # println(sum(σzτz[idx_chern_minus,:])/size(σzτz,2))
    println(length(idx_chern_minus)," ",length(idx_chern_plus))
    # maximal population of the Chern states, and random sprinkling of the remaining states (k unresolved)
    νmax = round(Int,(hf.ν+4)/8 * size(H0,1))
    states_to_populate = zeros(Int,νmax)
    # println(νmax," ",length(idx_chern_minus)+length(idx_chern_plus)/2)
    if νmax < length(idx_chern_minus)
        states_to_populate .= shuffle(idx_chern_minus)[1:νmax]
    else
        states_to_populate[1:length(idx_chern_minus)] .= idx_chern_minus
        states_to_populate[(length(idx_chern_minus)+1):end] .= shuffle(idx_chern_plus)[1:(νmax-length(idx_chern_minus))] 
        # id_running = length(idx_chern_minus)+1
        # for iplus in idx_chern_plus 
        #     if (iplus-1)÷(size(H0,1)÷4) +1 < 3
        #         # println(iplus)
        #         states_to_populate[id_running]=iplus
        #         id_running = id_running+1
        #     end
        # end
    end
    
    for ik in 1:size(hf.P,3)
        for iband in states_to_populate
            hf.P[:,:,ik] .+= conj(view(vecs,:,iband,ik))*transpose(view(vecs,:,iband,ik))
        end
    end

    for ik in 1:size(hf.P,3)
        hf.P[:,:,ik] .= hf.P[:,:,ik] - 0.5I
    end
    return nothing
end

function init_P_strong_coupling(hf::HartreeFock;
        P0::Array{ComplexF64,3}=ones(ComplexF64,1,1,1),H0::Array{ComplexF64,3}=ones(ComplexF64,1,1,1))
    # this function initializes the density matrix into a Chern state of the strong coupling spectrum 
    # first need to recreate the density matrix based on CNP 
    hf.H .= H0
    hf.P .= P0
    update_P(hf;α=1.0)
    println("Initialization based on populating excitation spectra of CNP")
    return nothing
end

function init_P_flavor_polarization(hf::HartreeFock)
    # maximal population of the Chern states, and random sprinkling of the remaining states 
    νmax = round(Int,(hf.ν+4)/8 * size(hf.H,1) * size(hf.H,3))
    n_per_valley_spin = size(hf.H,1)*size(hf.H,3)÷(hf.nη*hf.ns)

    tmpP = reshape(hf.P,hf.nb*hf.q,hf.nη*hf.ns,hf.nb*hf.q,hf.nη*hf.ns,:)

    num_full_flavors_to_populate = νmax ÷ n_per_valley_spin 
    num_part_flavors_to_populate = (νmax % n_per_valley_spin ==0) ? 0 : 1
    
    flavors_to_populate = randperm(hf.nη*hf.ns)[1:(num_full_flavors_to_populate+num_part_flavors_to_populate)]
    # flavors_to_populate = collect(1:hf.nη*hf.ns)[1:(num_full_flavors_to_populate+num_part_flavors_to_populate)]

    for iηs in flavors_to_populate[1:(end-num_part_flavors_to_populate)], iq in 1:size(tmpP,1)
        tmpP[iq,iηs,iq,iηs,:] .= 1.0 
    end
    if num_part_flavors_to_populate!=0
        iηs = flavors_to_populate[end]
        n_remain = νmax - (length(flavors_to_populate)-1)*n_per_valley_spin 
        states_to_populate = randperm(n_per_valley_spin)[1:n_remain]
        for ip in states_to_populate 
            iband = (ip-1)%(hf.nb*hf.q) + 1 
            ik = (ip-1)÷(hf.nb*hf.q) + 1
            tmpP[iband,iηs,iband,iηs,ik] = 1.0 
        end
    end
    
    for ik in 1:size(hf.P,3)
        hf.P[:,:,ik] .= hf.P[:,:,ik] - 0.5I
    end
    return nothing
end