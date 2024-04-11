include("BMLL_IKS.jl")

mutable struct QuantumGeometryBM
    # this contains information about the quantum geometric aspects 
    # of a magnetic subband given its wavefunctions 
    # only works for BM

    kvec::Array{ComplexF64}  # kvectors specifying the Brillouin zone 

    Uk::Array{ComplexF64,5}  # n x 2q x nq x q x nq , input to calculating geometric aspects 
    ϵk::Array{ComplexF64,3} # band basis, 2q x nq x nq 

    Q::Array{ComplexF64,5}  # α x β x length(kvec) x μ x ν  Non-abelian quantum geometric tensor 

    savename::String

    QuantumGeometryBM() = new()
end

function computeQuantumGeometryBM(bm::BMLL;savename::String="")
    qg = QuantumGeometryBM()

    qg.kvec = (reshape(bm.latt.k1[1:bm.nq],:,1,1)*bm.params.g1 .+
                reshape(bm.latt.k2[1:bm.nq],1,1,:)*bm.params.g2 .+ 
                reshape((0:(bm.q-1))./bm.q,1,:,1)*bm.params.g1)[:]
    
    qg.UK = permutedims(bm.vec,(1,2,4,3,5))
    qg.ϵk = bm.spectrum

    qg.Q = zeros(ComplexF64,size(bm.vec,2),size(bm.vec,2),length(qg.kvec),2,2)   # Q(m,n,k;x,y)

    computeBerryCurvature(qg)
    computeMetric(qg)

    if !isempty(savename)
        save(savename,"QuantumGeometryBM",qg)
    end

    return qg 
end

function computeBerryCurvature(qg::QuantumGeometryBM)
    return nothing 
end

function computeMetric(qg::QuantumGeometryBM)
    
    return nothing 
end