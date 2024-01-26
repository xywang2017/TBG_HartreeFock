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
q = 8 #parse(Int,ARGS[4])
ϕ = p//q
twist_angle = 1.20  # parse(Float64,ARGS[5])
_is_strain = "nostrain" # ARGS[6]

foldername =  @sprintf "NonInt/%d_%s" round(Int,twist_angle*100) _is_strain 
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
    println("p= ",p,", q= ",q,", nq= ",nq)
    fname = joinpath(fpath,"$(foldername)/_$(p)_$(q)_$(str)_metadata.jld2")
    println(fname)
    if isequal(_is_strain,"nostrain")
        params = Params(ϵ=0.00,Da=0.0,φ=0.0*π/180,dθ=twist_angle*π/180,w1=110,w0=110*w0,vf=2482)
    else
        params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*π/180,w1=110,w0=110*w0,vf=2482)
    end
    initParamsWithStrain(params)
    constructbmLL(bm,params;ϕ= ϕ,nLL=25*q÷p,nq=nq,fname=fname,α=w0, 
        _hBN=false,_strain=true, _σrotation=false, _valley=str,_calculate_overlap=true)

    return bm
end

#
bm = compute_bmLL(ϕ,str,w0,w0str);


function compute_mtg(bm::bmLL,ϕ::Rational,str::String,w0::Float64,w0str::String)
    p = numerator(ϕ)
    q = denominator(ϕ)
    if !isdir(joinpath(fpath,"$(foldername)"))
        mkpath(joinpath(fpath,"$(foldername)"))
    end
    fname = joinpath(fpath,"$(foldername)/_$(p)_$(q)_$(str)_mtg_metadata.jld2")
    mtg = constructMTG(bm;lr=10,fname=fname)
    return mtg
end

mtg = compute_mtg(bm,ϕ,str,w0,w0str);


ψr =  reshape(mtg.W[2,1,:],10,:)
rvec = reshape(mtg.coord.z,10,:)
fig = figure()
pcolormesh(real(rvec) ./abs(mtg.params.a1), imag(rvec)./abs(mtg.params.a1), angle.(ψr), cmap="hsv")
axis("equal")
tight_layout()
display(fig)
close(fig)

uvec = reshape(bm.vec,bm.nH*bm.p*2,2bm.q,bm.q*bm.nq^2);
ψ = zeros(ComplexF64,2,2bm.q,bm.q*bm.nq^2,mtg.coord.nr);
W = reshape(mtg.W,2,bm.nH*bm.p*2,bm.q*bm.nq^2,mtg.coord.nr);

for ik in 1:size(ψ,3), ib in 1:size(ψ,2), τ in 1:2 , ir in 1:mtg.coord.nr
    ψ[τ,ib,ir] = W[τ,:,ik,ir]'*uvec[:,ib,ik]
end

fig = figure()
pcolormesh(real(rvec) ./abs(mtg.params.a1), imag(rvec)./abs(mtg.params.a1), abs2.(reshape(ψ[1,1,1,:],10,:)), cmap="bwr")
axis("equal")
tight_layout()
display(fig)
close(fig)