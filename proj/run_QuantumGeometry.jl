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
# params = Params(ϵ=0.0,Da=0,φ=0.0*π/180,dθ=twist_angle*π/180,w1=110,w0=110*w0,vf=2482,δ=30.0)
initParamsWithStrain(params)
nq = 24÷q

# -------------------------- BM structure factor Related ----------------------- # 
if isequal(str,"K")
    fname = joinpath(fpath,"$(foldername)/_$(p)_$(q)_$(str)_metadata.jld2")
else
    fname = joinpath(fpath,"$(foldername)/_$(p)_$(q)_$(str)_$(q1)_$(q2)_metadata.jld2")
end
bm = bmLL();
constructbmLL(bm,params;ϕ=ϕ,nLL=25*q÷p,nq=nq,fname=fname,α=w0,
        _hBN=true,_strain=true, _σrotation=false, _valley=str,_calculate_overlap=true,q0=QIKS);
# -------------------------- Quantum Geometry Related ----------------------- # 
qg, tmpF, tmpG = computeQuantumGeometryBM(params;ϕ=ϕ,nq=nq,fname=fname,_valley=str,q0=QIKS);

# tmpF = readdlm("QuantumGeometry/_$(p)_$(q)_NonInt_F.txt")
# tmpG = readdlm("QuantumGeometry/_$(p)_$(q)_NonInt_G.txt")
# tmpF = [qg.F[ib,ib,ik] for ik in 1:size(qg.F,3) for ib in 1:size(qg.F,1)];
# tmpG = [qg.G[ib,ib,ik] for ik in 1:size(qg.F,3) for ib in 1:size(qg.F,1)];
# σF = sqrt(sum(abs2.(tmpF / (2π)*(q-p).- 1))/(nq^2*q*(q-p)) )
σF = sqrt(sum(abs2.(tmpF / (2π).- 1))/(nq^2*q) )
Tη = (sum(tmpG) - abs(sum(tmpF))) / (nq^2*q)
# ---------------------------- Berry curvature ----------------------- # 
fig,ax = subplots(2,1,figsize=(8,4))
lk = length(qg.latt.k1)
kvec = (reshape(qg.latt.k1,:,1) .+ 1im*reshape(qg.latt.k2,1,:))[:,1:qg.nq]
kvec = reshape(qg.latt.kvec,nq*q,:)[:,1:nq]/sqrt(abs(qg.params.g1)*abs(qg.params.g2))
pl= ax[1].pcolormesh(real(kvec),imag(kvec),reshape(tmpG,:,qg.nq),cmap="bwr")
colorbar(pl,shrink=0.7,ax=ax[1])
pl = ax[2].pcolormesh(real(kvec),imag(kvec),-reshape(tmpF,:,qg.nq),cmap="bwr")
colorbar(pl,shrink=0.7,ax=ax[2])
# for i in 1:2 
#     ax[i].axis("equal")
#     ax[i].axis("off")
#     ax[i].set_xticks([])
#     ax[i].set_yticks([])
#     for j in (1:(q-1))/q
#         ax[i].axvline(j-0.5/(nq*q),c="k",ls=":")
#     end
# end
tight_layout()
savefig("_$(p)_$(q)_strain_NonInt.png",transparent=true,dpi=500)
display(fig)
close(fig)

writedlm("_$(p)_$(q)_Ky.txt",imag(kvec))
# writedlm("_$(p)_$(q)_Int_G.txt",reshape(tmpG,:,qg.nq))
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