using PyPlot
using JLD2
fpath = joinpath(pwd(),"B0")
include(joinpath(fpath,"libs/HF_mod.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))

prefix =3
ν = 2.0
νstr = round(Int,1000*ν)
# ------------------ Specification ------------------ #
lk = 19
params = Params(ϵ=0.002,Da=-4100,dθ=1.06π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)
latt = Lattice()
initLattice(latt,params;lk=lk)

bm_path = joinpath(fpath,"data/strain_bm_lk$(lk).jld2")
hf_path = joinpath(fpath,"data/$(prefix)_strain_hf_$(νstr)_lk$(lk).jld2")
# hf_path = "typical_starting_point.jld2"

# ----------------- Hartree-Fock analysis part ---------------- # 
hf = load(hf_path,"hf");
kvec = reshape(latt.kvec ./ abs(params.g1),lk,lk)
ϵ0 = reshape(hf.ϵk,hf.nt,lk,lk)
plot_contour_maps(kvec,ϵ0[3,:,:],points=[params.Kt/abs(params.g1)],contourlines=[hf.μ])
iΓ = (lk%2==0) ? (lk÷2) : ((lk-1)÷2+1)
kcut = real(kvec[:,iΓ])
Ecut = ϵ0[:,:,iΓ]

plot_energy_cuts(kcut,Ecut,lines=[μ])


# ----------------- valley-spin polarization info ----------------- # 
s0 = ComplexF64[1 0;0 1]
s1 = ComplexF64[0 1;1 0]
s2 = ComplexF64[0 -1im;1im 0]
s3 = ComplexF64[1 0;0 -1]
ηz = zeros(Float64,size(hf.ϵk))
ηx = zeros(Float64,size(hf.ϵk))
ηy = zeros(Float64,size(hf.ϵk))
σz = zeros(Float64,size(hf.ϵk))
σx = zeros(Float64,size(hf.ϵk))
σy = zeros(Float64,size(hf.ϵk))
bz = zeros(Float64,size(hf.ϵk))
bx = zeros(Float64,size(hf.ϵk))
by = zeros(Float64,size(hf.ϵk))
c2T = zeros(Float64,size(hf.ϵk))
for ik in 1:size(hf.ϵk,2)
    F = eigen(Hermitian(view(hf.H,:,:,ik)))
    ηz[:,ik] = real(diag(F.vectors'*kron(s0,kron(s3,s0))*F.vectors))
    ηx[:,ik] = real(diag(F.vectors'*kron(s0,kron(s1,s0))*F.vectors))
    ηy[:,ik] = real(diag(F.vectors'*kron(s0,kron(s2,s0))*F.vectors))
    σz[:,ik] = real(diag(F.vectors'*kron(s3,kron(s0,s0))*F.vectors))
    σx[:,ik] = real(diag(F.vectors'*kron(s1,kron(s0,s0))*F.vectors))
    σy[:,ik] = real(diag(F.vectors'*kron(s2,kron(s0,s0))*F.vectors))
    bz[:,ik] = real(diag(F.vectors'*kron(s0,kron(s0,s3))*F.vectors))
    bx[:,ik] = real(diag(F.vectors'*kron(s0,kron(s0,s1))*F.vectors))
    by[:,ik] = real(diag(F.vectors'*kron(s0,kron(s0,s2))*F.vectors))
    c2T[:,ik] = real(diag(F.vectors'*conj.(F.vectors)))
end

# # C2T 
# plot_energy_cuts_with_order_parameters(kcut,Ecut,
#                 reshape(c2T,:,lk,lk)[:,:,iΓ],lines=[μ])

# #test valley coherence
plot_energy_cuts_with_order_parameters(kcut,Ecut,
                reshape(sqrt.(abs2.(ηx)+abs2.(ηy)),:,lk,lk)[:,:,iΓ],lines=[μ])

# # valley polarization
# plot_energy_cuts_with_order_parameters(kcut,Ecut,
#                 reshape((ηz.+1)/2,:,lk,lk)[:,:,iΓ],lines=[μ])

# # test spin polarization
plot_energy_cuts_with_order_parameters(kcut,Ecut,
                reshape((σz.+1)/2,:,lk,lk)[:,:,iΓ],lines=[μ])

# # test band mixing 
# plot_energy_cuts_with_order_parameters(kcut,Ecut,
#                 reshape(sqrt.(abs2.(bx)+abs2.(by)),:,lk,lk)[:,:,iΓ],lines=[μ])

strs, order_parameters = calculate_valley_spin_band_order_parameters(hf)
fig = figure(figsize=(2,10))
plot(order_parameters,eachindex(order_parameters),"b^")
yticks(collect(eachindex(order_parameters)),strs)
xlim([-1,2.5])
axvline(0,c="gray")
tight_layout()
display(fig)
close(fig)

# ----------------- Hartree-Fock statistics ---------------- # 
iter_oda = load(hf_path,"iter_oda");
iter_err = load(hf_path,"iter_err");
iter_energy = load(hf_path,"iter_energy");
println(iter_energy[end])
println(iter_err[end])

fig,ax = subplots(figsize=(4,3))
ax.plot(eachindex(iter_err),iter_err,"b-",label="err",ms=2)
ax1 = ax.twinx()
ax1.plot(eachindex(iter_energy),iter_energy,"r-",label="oda",ms=3)
ax.set_ylabel("Iter Err")
ax.set_yscale("log")
ax1.set_ylabel("Iter Energy")
# ax1.set_yscale("log")
tight_layout()
savefig("test.pdf")
display(fig)
close(fig)

