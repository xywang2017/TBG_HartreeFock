include("BMLL_IKS.jl")

mutable struct QuantumGeometryBM
    # this contains information about the quantum geometric aspects 
    # of a magnetic subband given its wavefunctions 
    # only works for BM

    fname::String  # BM structure factors contained here 

    _valley::String 

    params::Params
    latt::Lattice 

    p::Int 
    q::Int 
    nq::Int 

    δqs::Array{ComplexF64}  # q vectors to compute [Λ_q(k)]_α,β

    Λq::Array{ComplexF64,4}  # 2q x 2q x length(kvec) x n δqs , input to calculating geometric aspects
    U1::Array{ComplexF64,3} 
    
    F::Array{Float64,3}  # α x β x length(kvec)   Non-abelian quantum berry curvature
    G::Array{Float64,3}  # α x β x length(kvec)   Non-abelian quantum metric

    savename::String

    QuantumGeometryBM() = new()
end

function computeQuantumGeometryBM(params::Params;ϕ::Rational{Int}=1//10,U1::Array{ComplexF64,3}=zeros(ComplexF64,0,0,0),
                    nq::Int=2,fname::String="",savename::String="",_valley::String="K",q0::Complex{Int}=0+0im)
    qg = QuantumGeometryBM()

    qg.fname = fname 
    qg._valley = _valley
    qg.params = params 
    qg.p = numerator(ϕ)
    qg.q = denominator(ϕ)
    qg.nq = nq 
    qg.δqs = Complex{Int}[i+j*1im for i in -1:1 for j in -1:1]

    qg.latt = Lattice() 
    constructLatticeIKS(qg.latt,qg.params;lk = qg.nq*qg.q,_valley=qg._valley,q0=q0)

    qg.Λq = zeros(ComplexF64,2qg.q,2qg.q,qg.nq^2*qg.q,length(qg.δqs))
    if isempty(U1)
        qg.U1 = zeros(ComplexF64,2qg.q,2qg.q,qg.nq^2*qg.q)
        for ib in 1:(2qg.q) 
            qg.U1[ib,ib,:] .= 1.0+0.0im 
        end
    else
        qg.U1 = U1 
    end
    constructΛq(qg)

    tmpF = computeBerryCurvature(qg)
    tmpG = computeMetric(qg)

    if !isempty(savename)
        save(savename,"QuantumGeometryBM",qg)
    end

    return qg , tmpF, tmpG
end

function constructΛq(qg::QuantumGeometryBM)
    kvec = reshape( reshape(qg.latt.k1[1:qg.nq],:,1,1)*qg.params.g1 .+ 
            reshape(qg.latt.k2[1:qg.nq],1,1,:)*qg.params.g2 .+ 
            reshape((0:(qg.q-1))./qg.q,1,:,1)*qg.params.g1, qg.q*qg.nq,qg.nq )

    l1, l2 = size(kvec,1), size(kvec,2)

    Λ = zeros(ComplexF64,2qg.q*qg.nq^2*qg.q,2qg.q*qg.nq^2*qg.q)
    tmpΛ = reshape(Λ,2qg.q,qg.q,qg.nq,qg.nq,2qg.q,qg.q,qg.nq,qg.nq)
    tmpU1 = reshape(qg.U1,2qg.q,2qg.q,qg.q,qg.nq,qg.nq)
    for iq in eachindex(qg.δqs)
        δq1, δq2 = Int(real(qg.δqs[iq])), Int(imag(qg.δqs[iq]))
        for ik2 in 1:size(kvec,2), ik1 in 1:size(kvec,1)
            ik = (ik2-1)*l1 + ik1 
            ip1 = mod(ik1+δq1-1,l1) + 1
            ip2 = mod(ik2+δq2-1,l2) + 1
            m = ( (ik1+δq1) - ip1 ) ÷ l1 
            n = ( (ik2+δq2) - ip2 ) ÷ l2 
            jldopen(qg.fname,"r") do f 
                Λ = f["$(m)_$(n)"]
                _ik1, _rk = mod(ik1-1,qg.nq) + 1, div(ik1-1,qg.nq) + 1
                _ip1, _rp = mod(ip1-1,qg.nq) + 1, div(ip1-1,qg.nq) + 1
                tmpΛ = reshape(Λ,2qg.q,qg.q,qg.nq,qg.nq,2qg.q,qg.q,qg.nq,qg.nq)
                qg.Λq[:,:,ik,iq] = tmpU1[:,:,_rk,_ik1,ik2]'*
                                    tmpΛ[:,_rk,_ik1,ik2,:,_rp,_ip1,ip2]*tmpU1[:,:,_rp,_ip1,ip2]
            end
        end

    end
    
    return nothing 
end

function computeBerryCurvature(qg::QuantumGeometryBM)
    qg.F = zeros(Float64,2qg.q,2qg.q,qg.nq^2*qg.q)   # berry curvature
    δqs = [1+0im;0+1im;-1+0im;0-1im]
    iqs = [findfirst(x->x==δq,qg.δqs) for δq in δqs] 
    Λq = reshape(qg.Λq,2qg.q,2qg.q,qg.nq*qg.q,qg.nq,length(qg.δqs))

    idxF = qg.q
    tmpF = zeros(Float64,size(qg.Λq,3)) # single band

    δ = 1/(qg.nq*qg.q)
    for ik1 in 1:(qg.nq*qg.q), ik2 in 1:qg.nq
        ik = (ik2-1)*(qg.nq*qg.q) + ik1
        # F = Λq[:,:,ik1,ik2,iqs[1]]*
        #     Λq[:,:,mod(ik1+1-1,qg.nq*qg.q)+1,ik2,iqs[2]]*
        #     Λq[:,:,mod(ik1+1-1,qg.nq*qg.q)+1,mod(ik2+1-1,qg.nq)+1,iqs[3]]*
        #     Λq[:,:,ik1,mod(ik2+1-1,qg.nq)+1,iqs[4]]
        # qg.F[:,:,ik] = imag(log.(F)) ./δ^2 

        FF = Λq[idxF,idxF,ik1,ik2,iqs[1]]*
            Λq[idxF,idxF,mod(ik1+1-1,qg.nq*qg.q)+1,ik2,iqs[2]]*
            Λq[idxF,idxF,mod(ik1+1-1,qg.nq*qg.q)+1,mod(ik2+1-1,qg.nq)+1,iqs[3]]*
            Λq[idxF,idxF,ik1,mod(ik2+1-1,qg.nq)+1,iqs[4]]
        tmpF[ik] = - imag(log(FF))  /δ^2 
        # tmpF[ik] = imag(FF)  ./δ^2 
    end
    
    # above calculates δ^2 Im(⟨∂1u|∂2u⟩-⟨∂2u|∂1u⟩)
    # need to get Im(⟨∂xu|∂yu⟩-⟨∂yu|∂xu⟩)
    U = [real(qg.params.g1) real(qg.params.g2);imag(qg.params.g1) imag(qg.params.g2)]
    Uinv = inv(U)
    conversion_coeff = Uinv[1,1]*Uinv[2,2] - Uinv[2,1]*Uinv[1,2] 
    # conversion_coeff = 1
    qg.F .*= conversion_coeff * abs(imag(qg.params.g1'*qg.params.g2))/qg.q
    tmpF .*= conversion_coeff * abs(imag(qg.params.g1'*qg.params.g2))/qg.q
    return tmpF
end

function computeMetric(qg::QuantumGeometryBM)
    qg.G = zeros(Float64,2qg.q,2qg.q,qg.nq^2*qg.q)   # metric 
    δqs = [1+0im;-1+0im;0+1im;0-1im;1+1im;-1-1im]
    iqs = [findfirst(x->x==δq,qg.δqs) for δq in δqs] 
    Λq = reshape(qg.Λq,2qg.q,2qg.q,qg.nq*qg.q,qg.nq,length(qg.δqs))

    idxG = qg.q 
    tmpG = zeros(Float64,size(qg.Λq,3)) # single band

    U = [real(qg.params.g1) real(qg.params.g2);imag(qg.params.g1) imag(qg.params.g2)]
    Uinv = inv(U)
    Λ = Uinv*transpose(Uinv)
    δ = 1/(qg.nq*qg.q)
    for ik1 in 1:(qg.nq*qg.q), ik2 in 1:qg.nq
        ik = (ik2-1)*(qg.nq*qg.q) + ik1
        
        # g11 = ( I - Λq[:,:,ik1,ik2,iqs[1]] * 
        #                 Λq[:,:,mod(ik1+1-1,qg.nq*qg.q)+1,ik2,iqs[2]] ) ./ δ^2 
        # g22 = ( I - Λq[:,:,ik1,ik2,iqs[3]] * 
        #                 Λq[:,:,ik1,mod(ik2+1-1,qg.nq)+1,iqs[4]] ) ./ δ^2 
        # g12 = ( I - Λq[:,:,ik1,ik2,iqs[5]] * 
        #                 Λq[:,:,mod(ik1+1-1,qg.nq*qg.q)+1,mod(ik2+1-1,qg.nq)+1,iqs[6]] ) ./ δ^2 .- (g11+g22)
        
        # g12 ./= 2.0 

        # qg.G[:,:,ik] = Λ[1,1]*real(g11) + Λ[2,2]*real(g22) + Λ[1,2] * real(g12) * 2

        # g11 = ( I - Λq[idxG,idxG,ik1,ik2,iqs[1]] * 
        #                 Λq[idxG,idxG,mod(ik1+1-1,qg.nq*qg.q)+1,ik2,iqs[2]] ) / δ^2 
        # g22 = ( I - Λq[idxG,idxG,ik1,ik2,iqs[3]] * 
        #                 Λq[idxG,idxG,ik1,mod(ik2+1-1,qg.nq)+1,iqs[4]] ) / δ^2 
        # g12 = ( I - Λq[idxG,idxG,ik1,ik2,iqs[5]] * 
        #                 Λq[idxG,idxG,mod(ik1+1-1,qg.nq*qg.q)+1,mod(ik2+1-1,qg.nq)+1,iqs[6]] ) / δ^2 - (g11+g22)
        
        g11 = (I - Λq[idxG,idxG,ik1,ik2,iqs[1]] * 
                        Λq[idxG,idxG,mod(ik1+1-1,qg.nq*qg.q)+1,ik2,iqs[2]] ) / δ^2 
        g22 = (I - Λq[idxG,idxG,ik1,ik2,iqs[3]] * 
                        Λq[idxG,idxG,ik1,mod(ik2+1-1,qg.nq)+1,iqs[4]] ) / δ^2 
        g12 = (I - Λq[idxG,idxG,ik1,ik2,iqs[5]] * 
                        Λq[idxG,idxG,mod(ik1+1-1,qg.nq*qg.q)+1,mod(ik2+1-1,qg.nq)+1,iqs[6]] ) / δ^2 - (g11+g22)
        
        g12 /= 2.0 

        tmpG[ik] =  Λ[1,1]*real(g11) + Λ[2,2]*real(g22) + Λ[1,2] * real(g12) * 2
    end
    tmpG .*= abs(imag(qg.params.g1'*qg.params.g2))/qg.q
    return tmpG
end


function computeBerryCurvaturev2(qg::QuantumGeometryBM)
    # works for Dan Parker's result with hBN
    qg.F = zeros(Float64,qg.q-qg.p,qg.q-qg.p,qg.nq^2*qg.q)   # berry curvature
    δqs = [1+0im;0+1im;-1+0im;0-1im]
    iqs = [findfirst(x->x==δq,qg.δqs) for δq in δqs] 
    Λq = reshape(qg.Λq,2qg.q,2qg.q,qg.nq*qg.q,qg.nq,length(qg.δqs))

    idxF = qg.q 
    tmpF = zeros(Float64,size(qg.Λq,3)) # single band

    idxs = (qg.q+p+1):(2qg.q)
    δ = 1/(qg.nq*qg.q)
    for ik1 in 1:(qg.nq*qg.q), ik2 in 1:qg.nq
        ik = (ik2-1)*(qg.nq*qg.q) + ik1
        F = Λq[idxs,idxs,ik1,ik2,iqs[1]]*
            Λq[idxs,idxs,mod(ik1+1-1,qg.nq*qg.q)+1,ik2,iqs[2]]*
            Λq[idxs,idxs,mod(ik1+1-1,qg.nq*qg.q)+1,mod(ik2+1-1,qg.nq)+1,iqs[3]]*
            Λq[idxs,idxs,ik1,mod(ik2+1-1,qg.nq)+1,iqs[4]]
        tmp = svd(F)
        F = tmp.U*tmp.Vt 
        qg.F[:,:,ik] = -imag(log.(F)) ./δ^2 

        FF = Λq[idxF,idxF,ik1,ik2,iqs[1]]*
            Λq[idxF,idxF,mod(ik1+1-1,qg.nq*qg.q)+1,ik2,iqs[2]]*
            Λq[idxF,idxF,mod(ik1+1-1,qg.nq*qg.q)+1,mod(ik2+1-1,qg.nq)+1,iqs[3]]*
            Λq[idxF,idxF,ik1,mod(ik2+1-1,qg.nq)+1,iqs[4]]
        tmpF[ik] = - imag(log(FF))  /δ^2 
    end
    
    # above calculates δ^2 Im(⟨∂1u|∂2u⟩-⟨∂2u|∂1u⟩)
    # need to get Im(⟨∂xu|∂yu⟩-⟨∂yu|∂xu⟩)
    U = [real(qg.params.g1) real(qg.params.g2);imag(qg.params.g1) imag(qg.params.g2)]
    Uinv = inv(U)
    conversion_coeff = Uinv[1,1]*Uinv[2,2] - Uinv[2,1]*Uinv[1,2] 
    # conversion_coeff = 1
    qg.F .*= conversion_coeff 
    tmpF .*= conversion_coeff

    qg.F .*= abs(imag(qg.params.g1'*qg.params.g2))/qg.q
    return tmpF
end

function computeMetricv2(qg::QuantumGeometryBM)
    qg.G = zeros(Float64,qg.q-qg.p,qg.q-qg.p,qg.nq^2*qg.q)   # metric 
    δqs = [1+0im;-1+0im;0+1im;0-1im;1+1im;-1-1im]
    iqs = [findfirst(x->x==δq,qg.δqs) for δq in δqs] 
    Λq = reshape(qg.Λq,2qg.q,2qg.q,qg.nq*qg.q,qg.nq,length(qg.δqs))

    idxs = (qg.q+qg.p+1):(2qg.q)

    U = [real(qg.params.g1) real(qg.params.g2);imag(qg.params.g1) imag(qg.params.g2)]
    Uinv = inv(U)
    Λ = Uinv*transpose(Uinv)
    δ = 1/(qg.nq*qg.q)
    for ik1 in 1:(qg.nq*qg.q), ik2 in 1:qg.nq
        ik = (ik2-1)*(qg.nq*qg.q) + ik1
        
        g11 = (I-Λq[idxs,idxs,ik1,ik2,iqs[1]] * 
                        Λq[idxs,idxs,mod(ik1+1-1,qg.nq*qg.q)+1,ik2,iqs[2]] ) ./ δ^2 
        g22 = (I-Λq[idxs,idxs,ik1,ik2,iqs[3]] * 
                        Λq[idxs,idxs,ik1,mod(ik2+1-1,qg.nq)+1,iqs[4]] ) ./ δ^2 
        g12 = (I-Λq[idxs,idxs,ik1,ik2,iqs[5]] * 
                        Λq[idxs,idxs,mod(ik1+1-1,qg.nq*qg.q)+1,mod(ik2+1-1,qg.nq)+1,iqs[6]] ) ./ δ^2 .- (g11+g22)
        
        g12 ./= 2.0 

        qg.G[:,:,ik] = Λ[1,1]*real(g11) + Λ[2,2]*real(g22) + Λ[1,2] * real(g12) * 2
    end

    qg.G .*= abs(imag(qg.params.g1'*qg.params.g2))/q
    return nothing
end