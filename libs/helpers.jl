using LinearAlgebra
using Parameters
using DelimitedFiles
using GSL 
using ClassicalOrthogonalPolynomials

# contains model independent helper functions 
function inγ(iH::Int)
    if iH == 1
        n, iγ = 0, 0
    else
        n = (iH - 2) ÷ 2 + 1
        iγ = mod(iH - 2, 2) + 1
    end
    return n, 2iγ - 3
end

function projector_para(vec1::ComplexF64, vec2::ComplexF64)
    # projector of vec1 onto direction of vec2
    return real(conj(vec1) * vec2) / abs(vec2)
end

function projector_norm(vec1::ComplexF64, vec2::ComplexF64)
    # projector of vec1 onto norm of vec2,  e2xẑ
    return imag(conj(vec1) * vec2) / abs(vec2)
end

function _associatedlaguerre(n::Int, m::Int, cplus::ComplexF64, cminus::ComplexF64)
    x = -real(cplus * cminus)
    val = 0.0 + 0.0im
    if n >= m
        val = exp(-x / 2) * sqrt(factorial(big(m)) / factorial(big(n))) * big(cplus)^(n - m) * laguerrel(m, n - m, x)
    else
        val = exp(-x / 2) * sqrt(factorial(big(n)) / factorial(big(m))) * big(cminus)^(m - n) * laguerrel(n, m - n, x)
    end
    # if n >= m
    #     val = exp(-x / 2) * exp((sf_lnfact(m)-sf_lnfact(n))/2) * cplus^(n - m) * sf_laguerre_n(m, n - m, x)
    # else
    #     val = exp(-x / 2) * exp((sf_lnfact(n)-sf_lnfact(m))/2) * cminus^(m - n) * sf_laguerre_n(n, m - n, x)
    # end
    return (abs(val)<1e-16) ? ComplexF64(0.0) : ComplexF64(val)
end

function _associatedlaguerre_v1(nLL::Int, cplus::ComplexF64, cminus::ComplexF64)
    # first construct the matrix <n|c_-a + c_+ a^dagger |m>, diagonalize, then exponentiate
    # 0,1,.. nscale*nLL-1
    # x = -real(cplus*cminus)
    nscale = 1
    mat = zeros(ComplexF64, nscale*nLL, nscale*nLL)
    for i in 1:(size(mat,2)-1)
        mat[i+1, i] = sqrt(i) * cplus * (1im)
        mat[i, i+1] = sqrt(i) * cminus * (1im)
    end
    # if norm(mat - mat') > 1e-6
    #     println("error with hermitian")
    # end
    F = eigen(Hermitian(mat))
    # return view(F.vectors,1:nLL,:) * Diagonal(exp.(-1im * F.values)) * view(F.vectors,1:nLL,:)'
    return F.vectors * Diagonal(exp.(-1im * F.values)) * F.vectors'
end


function _associatedlaguerre_v2(nLL::Int, cplus::ComplexF64, cminus::ComplexF64)
    # first construct the matrix <n|c_-a + c_+ a^dagger |m>, diagonalize, then exponentiate
    # 0,1,.. nLL-1
    mat = zeros(ComplexF64, nLL, nLL)
    for n in 0:(nLL-1), m in 0:(nLL-1)
        mat[n+1,m+1] = _associatedlaguerre(n, m, cplus, cminus)
    end
    return mat
end


function _tLL_v1(T::Matrix{ComplexF64}, q::ComplexF64, nLL::Int, nH::Int,
    lB::Float64,θ0::Float64, θ1::Float64, θ2::Float64, _σrotation::Bool)
    # computes the matrix element of T exp(-iq r) in the LL basis, specified by (n1γ1l1,n2γ2l2)
    # T is 2x2 Hamiltonian in the sublattice basis
    # θ0 is a global rotation of axis in the presence of strain 

    cplus = -1im * lB / sqrt(2) * (real(q) - 1im * imag(q))
    cminus = -1im * lB / sqrt(2) * (real(q) + 1im * imag(q))
    # _alv1 = _associatedlaguerre_v1(nLL, cplus, cminus)
    _alv1 = _associatedlaguerre_v2(nLL, cplus, cminus)
    oLL = zeros(ComplexF64, nH, nH)
    for iH1 in 1:nH, iH2 in 1:nH
        n1, γ1 = inγ(iH1)
        n2, γ2 = inγ(iH2)
        if _σrotation == true
            if n1 != 0 && n2 != 0 # if not zeroth LL 
                oLL[iH1, iH2] = (T[1, 1] * (γ1 * γ2 * exp(1im * (θ2 - θ1))) * _alv1[n1, n2]
                                 + T[2, 2] * _alv1[n1+1, n2+1]
                                 + T[1, 2] * (1im * γ1 * exp(-1im * (θ1-θ0))) * _alv1[n1, n2+1]
                                 + T[2, 1] * (-1im * γ2 * exp(1im * (θ2-θ0))) * _alv1[n1+1, n2]) / 2.0
            elseif n1 == 0 && n2 == 0
                oLL[iH1, iH2] = T[2, 2] * _alv1[n1+1, n2+1]
            elseif n1 == 0 && n2 != 0
                oLL[iH1, iH2] = (T[2, 2] * _alv1[n1+1, n2+1]
                                 + T[2, 1] * (-1im * γ2 * exp(1im * (θ2-θ0))) * _alv1[n1+1, n2]) / sqrt(2.0)
            elseif n1 != 0 && n2 == 0
                oLL[iH1, iH2] = (T[2, 2] * _alv1[n1+1, n2+1]
                                 +T[1, 2] * (1im * γ1 * exp(-1im * (θ1-θ0))) * _alv1[n1, n2+1]) / sqrt(2.0)
            else
                println("Error with _tLL σrotation true")
            end
        else
            if n1 != 0 && n2 != 0 # if not zeroth LL 
                oLL[iH1, iH2] = (T[1, 1] * (γ1 * γ2) * _alv1[n1, n2]
                                 + T[2, 2] * _alv1[n1+1, n2+1]
                                 + T[1, 2] * (1im * γ1 * exp(1im*θ0)) * _alv1[n1, n2+1]
                                 + T[2, 1] * (-1im * γ2 * exp(-1im*θ0)) * _alv1[n1+1, n2]) / 2.0
            elseif n1 == 0 && n2 == 0
                oLL[iH1, iH2] = T[2, 2] * _alv1[n1+1, n2+1]
            elseif n1 == 0 && n2 != 0
                oLL[iH1, iH2] = (T[2, 2] * _alv1[n1+1, n2+1]
                                 + T[2, 1] * (-1im * γ2*exp(-1im*θ0)) * _alv1[n1+1, n2]) / sqrt(2.0)
            elseif n1 != 0 && n2 == 0
                oLL[iH1, iH2] = (T[2, 2] * _alv1[n1+1, n2+1]
                                 +  T[1, 2] * (1im * γ1*exp(1im*θ0)) * _alv1[n1, n2+1]) / sqrt(2.0)
            else
                println("Error with _tLL σrotation false")
            end
        end
    end
    # oLL *= exp(1im * real(q) * imag(q) * lB^2 / 2)
    return oLL
end


function _tLL(T::Matrix{ComplexF64}, q::ComplexF64, n1::Int, γ1::Int, n2::Int, γ2::Int,
    lB::Float64, θ1::Float64, θ2::Float64, _σrotation::Bool)
    # computes the matrix element of T exp(-iq r) in the LL basis, specified by (n1γ1l1,n2γ2l2)
    # T is 2x2 Hamiltonian in the sublattice basis 

    cplus = -1im * lB / sqrt(2) * (real(q) - 1im * imag(q))
    cminus = -1im * lB / sqrt(2) * (real(q) + 1im * imag(q))

    oLL = 0.0 + 0.0im
    if _σrotation == true
        if n1 != 0 && n2 != 0 # if not zeroth LL 
            oLL = (T[1, 1] * (γ1 * γ2 * exp(1im * (θ2 - θ1))) * _associatedlaguerre(n1 - 1, n2 - 1, cplus, cminus)
                   + T[2, 2] * _associatedlaguerre(n1, n2, cplus, cminus)
                   + T[1, 2] * (1im * γ1 * exp(-1im * θ1)) * _associatedlaguerre(n1 - 1, n2, cplus, cminus)
                   + T[2, 1] * (-1im * γ2 * exp(1im * θ2)) * _associatedlaguerre(n1, n2 - 1, cplus, cminus)) / 2.0
        elseif n1 == 0 && n2 == 0
            oLL = T[2, 2] * _associatedlaguerre(n1, n2, cplus, cminus)
        elseif n1 == 0 && n2 != 0
            oLL = (T[2, 2] * _associatedlaguerre(n1, n2, cplus, cminus)
                   +
                   T[2, 1] * (-1im * γ2 * exp(1im * θ2)) * _associatedlaguerre(n1, n2 - 1, cplus, cminus)) / sqrt(2.0)
        elseif n1 != 0 && n2 == 0
            oLL = (T[2, 2] * _associatedlaguerre(n1, n2, cplus, cminus)
                   +
                   T[1, 2] * (1im * γ1 * exp(-1im * θ1)) * _associatedlaguerre(n1 - 1, n2, cplus, cminus)) / sqrt(2.0)
        else
            println("Error with _tLL σrotation true")
        end
    else
        if n1 != 0 && n2 != 0 # if not zeroth LL 
            oLL = (T[1, 1] * γ1 * γ2 * _associatedlaguerre(n1 - 1, n2 - 1, cplus, cminus)
                   + T[2, 2] * _associatedlaguerre(n1, n2, cplus, cminus)
                   + T[1, 2] * (1im * γ1) * _associatedlaguerre(n1 - 1, n2, cplus, cminus)
                   + T[2, 1] * (-1im * γ2) * _associatedlaguerre(n1, n2 - 1, cplus, cminus)) / 2.0
        elseif n1 == 0 && n2 == 0
            oLL = T[2, 2] * _associatedlaguerre(n1, n2, cplus, cminus)
        elseif n1 == 0 && n2 != 0
            oLL = (T[2, 2] * _associatedlaguerre(n1, n2, cplus, cminus)
                   +
                   T[2, 1] * (-1im * γ2) * _associatedlaguerre(n1, n2 - 1, cplus, cminus)) / sqrt(2.0)
        elseif n1 != 0 && n2 == 0
            oLL = (T[2, 2] * _associatedlaguerre(n1, n2, cplus, cminus)
                   +
                   T[1, 2] * (1im * γ1) * _associatedlaguerre(n1 - 1, n2, cplus, cminus)) / sqrt(2.0)
        else
            println("Error with _tLL σrotation false")
        end
    end
    # oLL *= exp(1im * real(q) * imag(q) * lB^2 / 2)
    return oLL
end

function _tLLvalleyKprime(T::Matrix{ComplexF64}, q::ComplexF64, n1::Int, γ1::Int, n2::Int, γ2::Int,
    lB::Float64, θ1::Float64, θ2::Float64, _σrotation::Bool)
    # computes the matrix element of T exp(-iq r) in the LL basis, specified by (n1γ1l1,n2γ2l2)
    # T is 2x2 Hamiltonian in the sublattice basis 

    cplus = -1im * lB / sqrt(2) * (real(q) - 1im * imag(q))
    cminus = -1im * lB / sqrt(2) * (real(q) + 1im * imag(q))
    oLL = 0.0 + 0.0im
    if _σrotation == true
        if n1 != 0 && n2 != 0 # if not zeroth LL 
            oLL = (T[1, 1] * _associatedlaguerre(n1, n2, cplus, cminus)
                   + T[2, 2] * (γ1 * γ2 * exp(1im * (θ2 - θ1))) * _associatedlaguerre(n1 - 1, n2 - 1, cplus, cminus)
                   + T[1, 2] * (1im * γ2 * exp(1im * θ2)) * _associatedlaguerre(n1, n2 - 1, cplus, cminus)
                   + T[2, 1] * (-1im * γ1 * exp(-1im * θ1)) * _associatedlaguerre(n1 - 1, n2, cplus, cminus)) / 2.0
        elseif n1 == 0 && n2 == 0
            oLL = T[1, 1] * _associatedlaguerre(n1, n2, cplus, cminus)
        elseif n1 != 0 && n2 == 0
            oLL = (T[1, 1] * _associatedlaguerre(n1, n2, cplus, cminus)
                   +
                   T[2, 1] * (-1im * γ1 * exp(-1im * θ1)) * _associatedlaguerre(n1 - 1, n2, cplus, cminus)) / sqrt(2.0)
        elseif n1 == 0 && n2 != 0
            oLL = (T[1, 1] * _associatedlaguerre(n1, n2, cplus, cminus)
                   +
                   T[1, 2] * (1im * γ2 * exp(1im * θ2)) * _associatedlaguerre(n1, n2 - 1, cplus, cminus)) / sqrt(2.0)
        else
            println("Error with _tLL σrotation true")
        end
    else
        if n1 != 0 && n2 != 0 # if not zeroth LL 
            oLL = (T[1, 1] * _associatedlaguerre(n1, n2, cplus, cminus)
                   + T[2, 2] * (γ1 * γ2) * _associatedlaguerre(n1 - 1, n2 - 1, cplus, cminus)
                   + T[1, 2] * (1im * γ2) * _associatedlaguerre(n1, n2 - 1, cplus, cminus)
                   + T[2, 1] * (-1im * γ1) * _associatedlaguerre(n1 - 1, n2, cplus, cminus)) / 2.0
        elseif n1 == 0 && n2 == 0
            oLL = T[1, 1] * _associatedlaguerre(n1, n2, cplus, cminus)
        elseif n1 != 0 && n2 == 0
            oLL = (T[1, 1] * _associatedlaguerre(n1, n2, cplus, cminus)
                   +
                   T[2, 1] * (-1im * γ1) * _associatedlaguerre(n1 - 1, n2, cplus, cminus)) / sqrt(2.0)
        elseif n1 == 0 && n2 != 0
            oLL = (T[1, 1] * _associatedlaguerre(n1, n2, cplus, cminus)
                   +
                   T[1, 2] * (1im * γ2) * _associatedlaguerre(n1, n2 - 1, cplus, cminus)) / sqrt(2.0)
        else
            println("Error with _tLL σrotation true")
        end
    end
    # oLL *= exp(1im * real(q) * imag(q) * lB^2 / 2)
    return oLL
end

function _tLL_v1_valleyKprime(T::Matrix{ComplexF64}, q::ComplexF64, nLL::Int, nH::Int,
    lB::Float64, θ0::Float64,θ1::Float64, θ2::Float64, _σrotation::Bool)
    # computes the matrix element of T exp(-iq r) in the LL basis, specified by (n1γ1l1,n2γ2l2)
    # T is 2x2 Hamiltonian in the sublattice basis 
    # θ0 is a global rotation angle due to strain 
    
    cplus = -1im * lB / sqrt(2) * (real(q) - 1im * imag(q))
    cminus = -1im * lB / sqrt(2) * (real(q) + 1im * imag(q))
    # _alv1 = _associatedlaguerre_v1(nLL, cplus, cminus)
    _alv1 = _associatedlaguerre_v2(nLL, cplus, cminus)
    oLL = zeros(ComplexF64, nH, nH)
    for iH1 in 1:nH, iH2 in 1:nH
        n1, γ1 = inγ(iH1)
        n2, γ2 = inγ(iH2)
        if _σrotation == true
            if n1 != 0 && n2 != 0 # if not zeroth LL 
                oLL[iH1, iH2] = (T[1, 1] * _alv1[n1+1,n2+1]
                    + T[2, 2] * (γ1 * γ2 * exp(1im * (θ2 - θ1))) * _alv1[n1,n2]
                    + T[1, 2] * (1im * γ2 * exp(1im * (θ2-θ0))) * _alv1[n1+1,n2]
                    + T[2, 1] * (-1im * γ1 * exp(-1im * (θ1-θ0))) * _alv1[n1,n2+1]) / 2.0
            elseif n1 == 0 && n2 == 0
                oLL[iH1, iH2] = T[1, 1] * _alv1[n1+1,n2+1]
            elseif n1 != 0 && n2 == 0
                oLL[iH1, iH2] = (T[1, 1] * _alv1[n1+1,n2+1]
                    + T[2, 1] * (-1im * γ1 * exp(-1im * (θ1-θ0))) * _alv1[n1,n2+1]) / sqrt(2.0)
            elseif n1 == 0 && n2 != 0
                oLL[iH1, iH2] = (T[1, 1] * _alv1[n1+1,n2+1]
                    + T[1, 2] * (1im * γ2 * exp(1im * (θ2-θ0))) * _alv1[n1+1,n2]) / sqrt(2.0)
            else
                println("Error with _tLL σrotation true")
            end
        else
            if n1 != 0 && n2 != 0 # if not zeroth LL 
                oLL[iH1, iH2] = (T[1, 1] * _alv1[n1+1,n2+1]
                                + T[2, 2] * (γ1 * γ2 ) * _alv1[n1,n2]
                                + T[1, 2] * (1im * γ2*exp(-1im*θ0) ) * _alv1[n1+1,n2]
                                + T[2, 1] * (-1im * γ1*exp(1im*θ0)) * _alv1[n1,n2+1]) / 2.0
            elseif n1 == 0 && n2 == 0
                oLL[iH1, iH2] = T[1, 1] * _alv1[n1+1,n2+1]
            elseif n1 != 0 && n2 == 0
                oLL[iH1, iH2] = (T[1, 1] * _alv1[n1+1,n2+1]
                    + T[2, 1] * (-1im * γ1 * exp(1im*θ0) ) * _alv1[n1,n2+1]) / sqrt(2.0)
            elseif n1 == 0 && n2 != 0
                oLL[iH1, iH2] = (T[1, 1] * _alv1[n1+1,n2+1]
                    + T[1, 2] * (1im * γ2 *exp(-1im*θ0)) * _alv1[n1+1,n2]) / sqrt(2.0)
            else
                println("Error with _tLL σrotation false")
            end
        end
    end
    # oLL *= exp(1im * real(q) * imag(q) * lB^2 / 2)
    return oLL
end


function check_Hermitian(H::Matrix{ComplexF64})
    err = norm(H - H')
    if err > 1e-6
        println("Error with Hermitian Hamiltonian")
    end
    return nothing
end


function check_Unitary(A::Matrix{ComplexF64})
    err = norm(A'*A -I)
    if err > 1e-6
        println("Error with Unitarity of Matrix")
    end
    return nothing
end
