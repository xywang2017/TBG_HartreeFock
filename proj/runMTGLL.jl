using PyPlot
using Printf
using JLD2
fpath = pwd()
# include(joinpath(fpath,"libs/bmLL.jl"))
include(joinpath(fpath,"libs/mtg_real_space.jl"))

BLAS.set_num_threads(1)

str = "K" #ARGS[1]
w0 = 0.7 #parse(Float64,ARGS[2])*0.1
w0str = "07" #ARGS[2]
p = 1 #parse(Int,ARGS[3])
q = 6 #parse(Int,ARGS[4])
ϕ = p//q
twist_angle = 1.05  # parse(Float64,ARGS[5])
_is_strain = "strain" # ARGS[6]
q1 = 1
q2 = 0 
q0 = q1 + 1im*q2 

foldername =  @sprintf "MinHao/NonInt/%d_%s" round(Int,twist_angle*100) _is_strain 
# calculate spectrum
function compute_bmLL(ϕ::Rational,str::String,w0::Float64,w0str::String)
    p = numerator(ϕ)
    q = denominator(ϕ)
    #if !isdir(joinpath(fpath,"$(foldername)/B/data_w$(w0str)/_$(p)_$(q)"))
        #mkpath(joinpath(fpath,"$(foldername)/B/data_w$(w0str)/_$(p)_$(q)"))
    #end
    if !isdir(joinpath(fpath,"$(foldername)"))
        mkpath(joinpath(fpath,"$(foldername)"))
    end
    bm = bmLL()
    nq = (denominator(ϕ)>6) ? 1 : 2
    if q == 3 
        nq = 4 
    elseif q ==2 
        nq = 6
    end
    nq = 12÷q
    if q==7 
        nq = 2 
    end
    if nq ==0 
        nq = 1
    end
    println("p= ",p,", q= ",q,", nq= ",nq)
    fname = joinpath(fpath,"$(foldername)/_$(p)_$(q)_$(str)_metadata.jld2")
    fname = ""
    println(fname)
    if isequal(_is_strain,"nostrain")
        params = Params(ϵ=0.00,Da=0.0,φ=0.0*π/180,dθ=twist_angle*π/180,w1=110,w0=110*w0,vf=2482)
    else
        params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*π/180,w1=110,w0=110*w0,vf=2482)
    end
    initParamsWithStrain(params)
    constructbmLL(bm,params;ϕ= ϕ,nLL=25*q÷p,nq=nq,fname=fname,α=w0, 
        _hBN=false,_strain=true, _σrotation=false, _valley=str,_calculate_overlap=false)

    return bm
end

#
bm = compute_bmLL(ϕ,str,w0,w0str);


function compute_mtg(bm::bmLL,ϕ::Rational,str::String,w0::Float64,w0str::String,q0::Complex{Int}=0+0im)
    p = numerator(ϕ)
    q = denominator(ϕ)
    if !isdir(joinpath(fpath,"$(foldername)"))
        mkpath(joinpath(fpath,"$(foldername)"))
    end
    q1, q2 = real(q0), imag(q0)
    fname = joinpath(fpath,"$(foldername)/_$(p)_$(q)_mtg_$(q1)_$(q2)_metadata.jld2")
    mtg = constructMTG(bm;lr=16,fname=fname,q0=q0)
    return mtg
end

mtg = compute_mtg(bm,ϕ,str,w0,w0str,q0);


# rvec = reshape(mtg.coord.x1,:,1) .+ 1im*reshape(mtg.coord.x2,1,:)
rvec = reshape(mtg.coord.z,mtg.coord.lr,:)
uvec = reshape(bm.vec,bm.nH*bm.p,2,2bm.q,bm.q*bm.nq^2);
ψ = zeros(ComplexF64,2,2,2bm.q,bm.q*bm.nq^2,mtg.coord.nr);
W = reshape(mtg.W,2,bm.nH*bm.p,2,bm.q*bm.nq^2,mtg.coord.nr,2);

for ik in 1:size(ψ,4), ib in 1:size(ψ,3), τ in 1:2 , l in 1:2, ir in 1:mtg.coord.nr
    ψ[τ,l,ib,ik,ir] = transpose(W[τ,:,l,ik,ir,1])*uvec[:,l,ib,ik]
end

fig,ax = subplots(1,4,figsize=(8,4))
for i in 1:4 
    idx = 12
    vmax = maximum(abs2.(ψ[:,:,idx,1,:]))
    r,c = (i-1)%2 + 1, (i-1)÷2 + 1
    pl=ax[i].pcolormesh(real(rvec), imag(rvec), reshape(abs2.(ψ[r,c,idx,1,:]),mtg.coord.lr,:), cmap="bwr",vmin=0,vmax=vmax)
    ax[i].plot([0],[0],"k+")
    if i==4
        colorbar(pl,ax=ax[i],shrink=0.6)
    end
    ax[i].axis("equal")
end
tight_layout()
display(fig)
close(fig)


fig = figure(figsize=(8,8))
pl=pcolormesh(real(rvec), imag(rvec),reshape(sum(abs2.(ψ),dims=(1,2,3,4)),mtg.coord.lr,:), cmap="bwr")
plot([0],[0],"k+")
colorbar(pl,shrink=0.6)
axis("equal")
tight_layout()
display(fig)
close(fig)