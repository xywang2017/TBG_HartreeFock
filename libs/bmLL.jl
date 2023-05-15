include("helpers.jl")
include("params.jl")
include("latt.jl")
# using Arpack

mutable struct bmLL
    # this implements p/q sequence LL with strain

    # parameters 
    params::Params

    # flag for σmatrix rotation 
    _σrotation::Bool 
    # flag for hBN alignment 
    _hBN::Bool
    # flag for strain 
    _strain::Bool
    # String for valley 
    _valley::String

    p::Int 
    q::Int 
    qϕ::ComplexF64  # magnetic translation vector: ϕ/ϕ0 * 2π/L2 * ê_2
    nLL::Int   # 0,1,2... 
    nγ::Int    # 2, particle-hole symmetry 
    nH::Int    # nLL*nγ - 1

    lB::Float64  # magnetic length in absolute units 
    nq::Int # number of discrete momentum points in the range of [0,1/q)
    latt::Lattice

    qjs::Vector{Complex{Int}}
 
    H::Array{ComplexF64,8}  # diagonal + tunneling part 
    Σz::Array{ComplexF64,6}  #  σz operator 

    spectrum::Array{Float64,3}
    vec::Array{ComplexF64,5}
    PΣz::Array{ComplexF64,4}  #  Σz operator projected onto 2q states
    σz::Array{Float64,3} # σz in eigenbasis
    σzE::Array{Float64,3} # σz in energy eigenbasis

    fname::String  # file name to dump all calculations into

    ## Coulomb related 
    gvec::Array{ComplexF64,2}
    Λ::Array{ComplexF64,2}  # nk x mp x ig

    bmLL() = new()
end

function constructbmLL(A::bmLL,params::Params;
                        ϕ::Rational{Int}=1//10,nLL::Int=10,nq::Int=2,fname::String="placeholder.txt", α::Float64=0.3,
                        _σrotation::Bool = false, _hBN::Bool = false, _strain::Bool= true, _valley::String="K", _calculate_overlap::Bool = true)
    
    A._σrotation = _σrotation
    A._hBN = _hBN 
    A._strain = _strain
    A._valley = _valley
    A.p = numerator(ϕ)
    A.q = denominator(ϕ)
    A.nLL = nLL 
    A.nγ = 2
    A.nH = A.nLL * A.nγ - 1
    A.nq = nq
    
    A.params = params

    A.qϕ = 2π/abs2(A.params.a2) * A.params.a2 * A.p/A.q
    
    A.latt = Lattice() 
    constructLattice(A.latt,A.params;lk = A.nq*A.q)  #[0,1)x[0,1), so far works for p/q < 1

    A.lB = sqrt( A.q/(2π*abs(A.p)) * A.params.area )
    A.qjs = (isequal(A._valley,"K")) ? Complex{Int}[-1; -1+1im; 1im] : Complex{Int}[1; 1-1im; -1im]
    # A.qjs = (isequal(A._valley,"K")) ? Complex{Int}[0; 1im; 1+1im] : Complex{Int}[0; -1im; -1-1im]
    
    A.H = zeros(ComplexF64,A.nH,A.p,2,A.nH,A.p,2,A.nq,A.nq)
    constructDiagonals(A)
    constructOffDiagonals(A)
    constructΣz(A)
    computeSpectrum(A)
    enforceC2P(A)
    # computeSpectrum_remote(A)

    A.fname = fname 
    # write energies and wavefunctions
    jldopen(fname, "w") do file
        file["E"] = A.spectrum 
        file["Vec"] = A.vec
        file["PΣz"] = A.PΣz
    end

    ng = 3
    A.gvec = reshape(collect(-ng:ng),:,1)*A.params.g1 .+ reshape(collect(-ng*A.q:ng*A.q),1,:)*A.params.g2 ./A.q
    A.Λ = zeros(ComplexF64,2A.q*A.q*A.nq^2,2A.q*A.q*A.nq^2)
    
    if _calculate_overlap
        # for m in -ng:ng, n in -ng*A.q:ng*A.q 
        for m in -ng:ng, n in (ng*A.q):-1:(-ng*A.q)
            println("m:",m," n:",n)
            computeCoulombOverlap(A,m,n)
            jldopen(fname, "a") do file
                file["$(m)_$(n)"] = A.Λ
            end
        end
    end

    return nothing
end

function constructDiagonals(A::bmLL)
    ϵB = A.params.vf / A.lB
    for iH in 1:A.nH
        n,γ = inγ(iH)
        for ip in 1:A.p
            A.H[iH,ip,1,iH,ip,1,:,:] .-= γ*sqrt(2n) * ϵB
            A.H[iH,ip,2,iH,ip,2,:,:] .-= γ*sqrt(2n) * ϵB
        end
    end

    if (A._hBN==true)
        sign_hBN = isequal(A._valley,"K") ? 1 : (-1)
        for ip in 1:A.p
            A.H[1,ip,1,1,ip,1,:,:] .+= -A.params.δ * sign_hBN
            for n in 1:(A.nLL-1)
                for iγ1 in 1:2
                    γ1 = 2iγ1 - 3
                    iH1 = (n-1)*2 + iγ1 + 1 
                    for iγ2 in 1:2
                        γ2 = 2iγ2 - 3
                        iH2 = (n-1)*2 + iγ2 + 1 
                        A.H[iH1,ip,1,iH2,ip,1,:,:] .+= -(1-γ1*γ2)/2 * A.params.δ * sign_hBN
                    end
                end
            end
        end
    end

    if (A._strain==true)
        for iH in 1:A.nH
            n,γ = inγ(iH)
            def_pot = A.params.Da * (A.params.S[1,1]+A.params.S[2,2])/2
            for ip in 1:A.p
                A.H[iH,ip,1,iH,ip,1,:,:] .-= def_pot
                A.H[iH,ip,2,iH,ip,2,:,:] .+= def_pot
            end
        end
    end

    return nothing
end

function constructOffDiagonals(A::bmLL)
    r2 = 0
    _q = 0.0 + 0.0im
    θ_strain = angle(A.params.a2) - π/2
    T12 = zeros(ComplexF64,A.nH,A.nH)
    for j in eachindex(A.qjs)
        if (j==1)
            _Tj = A.params.T0 
        elseif (j==2)
            _Tj = A.params.T1
        else 
            _Tj = A.params.T2
        end
        _q = real(A.qjs[j])*A.params.g1 + imag(A.qjs[j])*A.params.g2
        _q = isequal(A._valley,"K") ? (_q + A.params.Kb - A.params.Kt) : (_q + A.params.Kt - A.params.Kb)
        # rewrite it as qx + i qy in rotated coordinate system, important for strained calculation
        _q = projector_norm(_q,A.params.a2) + 1im * projector_para(_q,A.params.a2)  
        if isequal(A._valley,"K")
            T12 .= _tLL_v1(_Tj,_q,A.nLL,A.nH,A.lB,-A.params.dθ/2,A.params.dθ/2,A._σrotation)
        else
            T12 .= _tLL_v1_valleyKprime(conj(_Tj),_q,A.nLL,A.nH,A.lB,-A.params.dθ/2,A.params.dθ/2,A._σrotation)
        end
        # for iH1 in 1:A.nH, iH2 in 1:A.nH
        #     n1,γ1 = inγ(iH1)
        #     n2,γ2 = inγ(iH2)
        #     T12[iH1,iH2] = _tLL(_Tj,_q,n1,γ1,n2,γ2,A.lB,-A.params.dθ/2,A.params.dθ/2,A._σrotation)
        #     # T12[iH1,iH2] = _tLL(_Tj,_q,n1,γ1,n2,γ2,A.lB,θ_strain,θ_strain,true)
        # end
        Kl = isequal(A._valley,"K") ? A.params.Kb : A.params.Kt
        Kr = isequal(A._valley,"K") ? A.params.Kt : A.params.Kb
        for ik2 in 1:A.nq, r1 in 0:(A.p-1)
            k2l = projector_para(A.params.g2,A.params.a2) * A.latt.k2[r1*A.nq+ik2] - projector_para(Kl,A.params.a2)
            r2 = mod(r1 + A.q * imag(A.qjs[j]),A.p)
            s = - ((r1 + A.q * imag(A.qjs[j])) - r2 )÷A.p
            p2 = A.latt.k2[r2*A.nq+ik2]
            for ik1 in 1:A.nq
                expfactor = exp(1im * 2π * s * (A.latt.k1[ik1]-p2*projector_para(A.params.a1,A.params.a2)/abs(A.params.a2)) ) * 
                            exp(1im *s*(s-1)/2 *projector_para(A.qϕ,A.params.a1)*abs(A.params.a1)) * 
                            exp(-1im * s * projector_norm(Kr,A.params.a2)*projector_norm(A.params.a1,A.params.a2)) * 
                            exp(1im * real(_q)*k2l*A.lB^2) * exp(1im * real(_q) * imag(_q) * A.lB^2 / 2)
                A.H[:,r1+1,1,:,r2+1,2,ik1,ik2] .+=  expfactor * T12
            end
        end
    end
    return nothing
end

function constructΣz(A::bmLL)
    A.Σz = zeros(ComplexF64,A.nH,A.p,2,A.nH,A.p,2)
    sign_Σz = isequal(A._valley,"K") ? 1 : -1
    for ip in 1:A.p
        A.Σz[1,ip,1,1,ip,1] = -1 * sign_Σz
        A.Σz[1,ip,2,1,ip,2] = -1 * sign_Σz
        for n in 1:(A.nLL-1)
            for iγ1 in 1:2
                γ1 = 2iγ1 - 3
                iH1 = (n-1)*2 + iγ1 + 1 
                for iγ2 in 1:2
                    γ2 = 2iγ2 - 3
                    iH2 = (n-1)*2 + iγ2 + 1 
                    A.Σz[iH1,ip,1,iH2,ip,1] = -(1-γ1*γ2)/2 * sign_Σz
                    A.Σz[iH1,ip,2,iH2,ip,2] = -(1-γ1*γ2)/2 * sign_Σz
                end
            end
        end
    end
    
    return nothing
end

function computeSpectrum(A::bmLL)
    A.vec = zeros(ComplexF64,2A.nH*A.p,2A.q,A.q,A.nq,A.nq) # inner x 2q x rk1 x nkq x nk1
    A.spectrum = zeros(Float64,2A.q,A.nq,A.nq)
    A.PΣz = zeros(ComplexF64,2A.q,2A.q,A.nq,A.nq)
    A.σz = zeros(Float64,2A.q,A.nq,A.nq)
    A.σzE = zeros(Float64,2A.q,A.nq,A.nq)

    Σz = reshape(A.Σz,2A.nH*A.p,2A.nH*A.p)
    H = reshape(A.H,2A.nH*A.p,2A.nH*A.p,A.nq,A.nq)
    for i2 in 1:size(H,4), i1 in 1:size(H,3)
        idmid = A.nH*A.p
        idx_flat = (idmid+1-A.q):(idmid+A.q)
        # F = eigen( Hermitian( H[:,:,i1,i2], :U))
        # A.spectrum[:,i1,i2] = F.values[idx_flat]
        # if norm(F.vectors*F.vectors' - I) >1e-8
        #     println("Error with normalization LL")
        # end  
        # # pick out 2q states in the middle 
        # vec .= F.vectors[:,idx_flat]
        
        # vals,vec = eigs(Hermitian(view(H,:,:,i1,i2), :U),nev=2A.q,which=:SM)
        vals,vec = eigen(Hermitian(view(H,:,:,i1,i2), :U),idx_flat)
        # if norm(imag(vals))>1e-6
        #     println("Error with Hermitican",norm(imag(vals))," ",i1," ",i2)
        # end
        idx_perm = sortperm(real(vals))
        A.spectrum[:,i1,i2]= real(vals[idx_perm])
        A.vec[:,:,1,i1,i2] = @view vec[:,idx_perm]
        generateU(A.vec,A.q,A.p,A.nq)
        A.PΣz[:,:,i1,i2] .= view(vec,:,idx_perm)' * Σz * view(vec,:,idx_perm)
        if imag(tr(A.PΣz[:,:,i1,i2]))>1e-6
            println("Error with realness of tr(PΣz)")
        end
        A.σzE[:,i1,i2] = real(diag(A.PΣz[:,:,i1,i2]))[idx_perm]
        A.σz[:,i1,i2] = real(eigvals(Hermitian(A.PΣz[:,:,i1,i2])))
    end
    return nothing
end

function enforceC2P(A::bmLL)
    return nothing
end

function computeSpectrum_remote(A::bmLL)
    A.vec = zeros(ComplexF64,2A.nH*A.p,4A.q,A.q,A.nq,A.nq) # inner x 2q x rk1 x nkq x nk1
    A.spectrum = zeros(Float64,4A.q,A.nq,A.nq)
    A.PΣz = zeros(ComplexF64,4A.q,4A.q,A.nq,A.nq)
    A.σz = zeros(Float64,4A.q,A.nq,A.nq)
    A.σzE = zeros(Float64,4A.q,A.nq,A.nq)

    Σz = reshape(A.Σz,2A.nH*A.p,2A.nH*A.p)
    H = reshape(A.H,2A.nH*A.p,2A.nH*A.p,A.nq,A.nq)
    for i2 in 1:size(H,4), i1 in 1:size(H,3)
        idmid = A.nH*A.p
        idx_flat = (idmid+1-2A.q):(idmid+2A.q)
        # F = eigen( Hermitian( H[:,:,i1,i2], :U))
        # A.spectrum[:,i1,i2] = F.values[idx_flat]
        # if norm(F.vectors*F.vectors' - I) >1e-8
        #     println("Error with normalization LL")
        # end  
        # # pick out 2q states in the middle 
        # vec .= F.vectors[:,idx_flat]
        
        # vals,vec = eigs(Hermitian(view(H,:,:,i1,i2), :U),nev=2A.q,which=:SM)
        vals,vec = eigen(Hermitian(view(H,:,:,i1,i2), :U),idx_flat)
        # if norm(imag(vals))>1e-6
        #     println("Error with Hermitican",norm(imag(vals))," ",i1," ",i2)
        # end
        idx_perm = sortperm(real(vals))
        A.spectrum[:,i1,i2]= real(vals[idx_perm])
        A.vec[:,:,1,i1,i2] = @view vec[:,idx_perm]
        # generateU(A.vec,A.q,A.p,A.nq)
        A.PΣz[:,:,i1,i2] .= view(vec,:,idx_perm)' * Σz * view(vec,:,idx_perm)
        if imag(tr(A.PΣz[:,:,i1,i2]))>1e-6
            println("Error with realness of tr(PΣz)")
        end
        A.σzE[:,i1,i2] = real(diag(A.PΣz[:,:,i1,i2]))[idx_perm]
        A.σz[:,i1,i2] = real(eigvals(Hermitian(A.PΣz[:,:,i1,i2])))
    end
    return nothing
end


function data_writeout(A::bmLL)
    writedlm(A.fname,[A.spectrum[:] A.σzE[:] A.σz[:]])
    return nothing
end

function generateU(vec::Array{ComplexF64,5},q::Int,p::Int,nq::Int)
    #This function generates narrow band eigenstates 
    uvec = reshape(view(vec,:,:,1,:,:),:,p,2,2q,nq,nq)
    rk2 = collect(0:(p-1))
    for r1 in 1:(q-1)
        i1 = 0
        while i1<q
            if mod((i1*p),q) == r1
                break 
            else 
                i1 = i1+1 
            end
        end
        vec[:,:,r1+1,:,:] .= reshape(uvec .* reshape(exp.(-1im*2π*i1*(rk2./q) ),1,:,1,1,1,1),:,2q,nq,nq)
    end
    return nothing
end

function computeCoulombOverlap(A::bmLL,m::Int,n::Int)
    # q = p - k + m g1 + n/q * g2, m from -ng to ng, n from -ng*q to ng*q
    # p,k in range [0,1)x[0,1/q), and on a nq x nq grid 
    A.Λ .= 0.0 + 0.0im
    tmpΛ = reshape(A.Λ,2A.q,A.q,A.nq,A.nq,2A.q,A.q,A.nq,A.nq)
    ΛΨ = zeros(ComplexF64,2A.nH*A.p,2A.nH*A.p)
    tmpΛΨ = reshape(ΛΨ,A.nH,A.p,2,A.nH,A.p,2)
    Ov = ComplexF64[1 0;0 1]
    for ip2 in 1:A.nq, ik2 in 1:A.nq 
        for ip1 in 1:A.nq, rp1 in 0:(A.q-1), ik1 in 1:A.nq, rk1 in 0:(A.q-1)
            _k1, _p1 = A.latt.k1[ik1+rk1*A.nq], A.latt.k1[ip1+rp1*A.nq]
            for rk2 in 0:(A.p-1)
                rp2, s = mod(rk2+n,A.p), -((rk2+n)-mod(rk2+n,A.p))÷A.p
                _k2, _p2 = A.latt.k2[ik2+rk2*A.nq], A.latt.k2[ip2+rp2*A.nq]
                _k20, _p20 = A.latt.k2[ik2], A.latt.k2[ip2]
                _q = (_p1-_k1 + m)*A.params.g1 +  (_p20 - _k20 + n/A.q)*A.params.g2
                _q = projector_norm(_q,A.params.a2) + 1im * projector_para(_q,A.params.a2)  

                for layer in 1:2
                    if isequal(A._valley,"K")
                        Kl = (layer==1) ? A.params.Kb : A.params.Kt
                    else
                        Kl = (layer==1) ? A.params.Kt : A.params.Kb
                    end
                    θl = (2layer-3)*A.params.dθ/2
                    k2l = projector_para(A.params.g2,A.params.a2) * _k2 - projector_para(Kl,A.params.a2)
                    expfactor = exp(1im * 2π * s * (_p1-_p2*projector_para(A.params.a1,A.params.a2)/abs(A.params.a2)) ) * 
                                    exp(1im *s*(s-1)/2 *projector_para(A.qϕ,A.params.a1)*abs(A.params.a1)) * 
                                    exp(-1im * s * projector_norm(Kl,A.params.a2)*projector_norm(A.params.a1,A.params.a2)) * 
                                    exp(1im * real(_q)*k2l*A.lB^2) * exp(1im * real(_q) * imag(_q) * A.lB^2 / 2)

                    if isequal(A._valley,"K")
                        tmpΛΨ[:,rk2+1,layer,:,rp2+1,layer] .= _tLL_v1(Ov,_q,A.nLL,A.nH,A.lB,θl,θl,A._σrotation) * expfactor
                    else
                        tmpΛΨ[:,rk2+1,layer,:,rp2+1,layer] .= _tLL_v1_valleyKprime(Ov,_q,A.nLL,A.nH,A.lB,θl,θl,A._σrotation) * expfactor
                    end
                end
            end
            tmpΛ[:,rk1+1,ik1,ik2,:,rp1+1,ip1,ip2] = view(A.vec,:,:,rk1+1,ik1,ik2)' * ΛΨ * view(A.vec,:,:,rp1+1,ip1,ip2)
        end
    end
    return nothing
end
