using PyPlot
using Printf
using JLD2
fpath = pwd()
include(joinpath(fpath,"libs/QuantumGeometry.jl"))

BLAS.set_num_threads(1)

str = "K"
w0, w0str = 0.7, "07"
p, q = 1, 4
ϕ = p//q
twist_angle = 1.05
_is_strain = "strain"
q1, q2  = 0, 0 
QIKS = q1 + 1im*q2 

foldername =  @sprintf "NonInt/%d_%s" round(Int,twist_angle*100) _is_strain 
if !isdir(joinpath(fpath,"$(foldername)"))
    mkpath(joinpath(fpath,"$(foldername)"))
end
if isequal(_is_strain,"nostrain")
    params = Params(ϵ=0.00,Da=0.0,φ=0.0*π/180,dθ=twist_angle*π/180,w1=110,w0=110*w0,vf=2482)
else
    params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*π/180,w1=110,w0=110*w0,vf=2482)
end
initParamsWithStrain(params)
nq = 6 #12÷q

# -------------------------- BM structure factor Related ----------------------- # 
if isequal(str,"K")
    fname = joinpath(fpath,"$(foldername)/_$(p)_$(q)_$(str)_metadata.jld2")
else
    fname = joinpath(fpath,"$(foldername)/_$(p)_$(q)_$(str)_$(q1)_$(q2)_metadata.jld2")
end
bm = bmLL();
constructbmLL(bm,params;ϕ=ϕ,nLL=25*q÷p,nq=nq,fname=fname,α=w0, 
        _hBN=false,_strain=true, _σrotation=false, _valley=str,_calculate_overlap=true,q0=QIKS);
# -------------------------- Quantum Geometry Related ----------------------- # 
qg = computeQuantumGeometryBM(params;ϕ=ϕ,nq=nq,fname=fname,_valley=str,q0=QIKS);



# ---------------------------- Berry curvature ----------------------- # 
fig = figure(figsize=(4,3))
kvec = reshape( reshape(qg.latt.k1[1:qg.nq],:,1,1)*qg.params.g1 .+ 
    reshape(qg.latt.k2[1:qg.nq],1,1,:)*qg.params.g2 .+ 
    reshape((0:(qg.q-1))./qg.q,1,:,1)*qg.params.g1, qg.q*qg.nq,qg.nq ) ./abs(params.g1)

# pcolormesh(imag(kvec),real(kvec),reshape(qg.F[q,q,:],:,qg.nq)/π,cmap="bwr")
imshow(reshape((qg.F[1,1,:]),:,qg.nq),origin="lower",extent=(1,nq+1,1,q*nq+1).-0.5)
# axis("equal")
colorbar()
tight_layout()
display(fig)
close(fig)

[sum(qg.F[i,i,:]) for i in 1:(2q)]
# ---------------------------- Structure factor ----------------------- # 
fig = figure(figsize=(8,2))
imshow(reshape(abs.(qg.Λq[qg.q,qg.q,:,9]),:,qg.nq)',origin="lower",extent=(1,nq+1,1,q*nq+1).-0.5,cmap="bwr")
for i in [4,8,12]
    axvline(i+0.5,ls=":",c="k")
end
colorbar(shrink=0.5)
tight_layout()
display(fig)
close(fig)