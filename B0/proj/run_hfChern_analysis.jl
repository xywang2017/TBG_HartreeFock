using PyPlot
using JLD2
fpath = pwd() 
include(joinpath(fpath,"B0/libs/HFChern_mod.jl"))
include(joinpath(fpath,"B0/libs/plot_helpers.jl"))

prefix = 1
flag = "random"
# νs = collect(0.0:0.2:4.0)
ν = 2.0
νstr = round(Int,1000*ν)
# ------------------ Specification ------------------ #
lk = 17
# params = Params(ϵ=0.00,Da=0,dθ=1.06π/180,w1=110,w0=77,vf=2482)
params = Params(ϵ=0.002,Da=-4100,dθ=1.05π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)
latt = Lattice()
initLattice(latt,params;lk=lk)

bm_path = joinpath(fpath,"feldman/B0/data/bm_lk$(lk).jld2")
hf_path = joinpath(fpath,"feldman/B0/data/$(prefix)_$(flag)_hf_$(νstr)_lk$(lk).jld2")
# ----------------- Hartree-Fock analysis part ---------------- # 
hf = load(hf_path,"hf");
iter_energy = load(hf_path,"iter_energy");
println(iter_energy[end])
kvec = reshape(latt.kvec ./ abs(params.g1),lk,lk)
ϵ0 = reshape(hf.ϵk,hf.nt,lk,lk)
# plot_contour_maps(kvec,ϵ0[5,:,:],points=[params.Kt/abs(params.g1)],contourlines=[hf.μ])
iΓ = (lk%2==0) ? (lk÷2) : ((lk-1)÷2+1)
kcut = real(kvec[:,iΓ])
Ecut = ϵ0[:,:,iΓ]
# plot_energy_cuts(kcut,Ecut,lines=[hf.μ])

# ----------------- valley-spin-bamd polarization info ----------------- # 
fig, ax = subplots(1,2,sharex=true,figsize=(4,5))
s = length(hf.Δ)
ax[1].axvline(0,c="gray")
ax[1].plot(hf.Δ[2:s÷2],eachindex(hf.Δ)[2:s÷2],"b^")
ax[1].set_yticks(collect(eachindex(hf.Δ))[2:s÷2])
ax[1].set_yticklabels(hf.Δstr[2:s÷2])
ax[2].axvline(0,c="gray")
ax[2].plot(hf.Δ[(s÷2+1):end],eachindex(hf.Δ)[(s÷2+1):end],"b^")
ax[2].set_yticks(collect(eachindex(hf.Δ))[(s÷2+1):end])
ax[2].set_yticklabels(hf.Δstr[(s÷2+1):end])
ax[2].set_xlim(-0.1,0.1)
tight_layout()
savefig("test.pdf")
display(fig)
close(fig)

s0 = ComplexF64[1 0;0 1]
s1 = ComplexF64[0 1;1 0]
s2 = ComplexF64[0 -1im;1im 0]
s3 = ComplexF64[1 0;0 -1]
Δ = zeros(size(hf.ϵk))
for ik in 1:size(hf.ϵk,2)
    F = eigen(Hermitian(view(hf.H,:,:,ik)))
    Δ[:,ik] = real(diag(F.vectors'*kron(s1,kron(s0,s3))*F.vectors))
end
# Δ .= hf.σzτz

plot_energy_cuts_with_order_parameters(kcut,Ecut,
                reshape(Δ,:,lk,lk)[:,:,iΓ],lines=[hf.μ])


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

### all the chemical potentials 
νs = collect(0:0.2:4.0)
μs = Float64[]
actual_νs = Float64[]
for ν in νs 
    νstr = round(Int,1000*ν)
    hf_path = joinpath(fpath,"data/2_strain_hf_$(νstr)_lk19.jld2")
    if ispath(hf_path)
        hf = load(hf_path,"hf");
        push!(actual_νs,round(Int,(hf.ν+4)/8*size(hf.H,1)*size(hf.H,3))/(size(hf.H,1)*size(hf.H,3))*8 - 4)
        push!(μs,hf.μ)
    end
end
fig = figure(figsize=(5,3))
plot(sort(actual_νs),μs[sortperm(actual_νs)],"b-o",ms=3)
xlabel("ν")
ylabel("μ")
axhline(0,ls="--",c="gray")
# ylim([0,30])
# ylim([0,14])
xlim([-0.1,4.1])
tight_layout()
# savefig("cascade_strain_gapless.pdf")
display(fig)
close(fig)

