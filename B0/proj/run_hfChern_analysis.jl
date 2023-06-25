using PyPlot
using JLD2
fpath = pwd() 
include(joinpath(fpath,"B0/libs/HFChern_mod.jl"))
include(joinpath(fpath,"B0/libs/plot_helpers.jl"))

prefix = 1
flag = "sp"
phi = 0
strain = 0
# νs = collect(0.0:0.2:4.0)
ν = -2.0
νstr = round(Int,1000*ν)
# ------------------ Specification ------------------ #
lk = 20
# params = Params(ϵ=0.00,Da=0,dθ=1.06π/180,w1=110,w0=77,vf=2482)
params = Params(ϵ=0.001*strain,Da=-4100,dθ=1.05π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)
latt = Lattice()
initLattice(latt,params;lk=lk)

bm_path = joinpath(fpath,"princeton/B0/data/strain$(strain)/phi$(phi)/bm_lk$(lk).jld2")
hf_path = joinpath(fpath,"princeton/B0/data/strain$(strain)/phi$(phi)/$(prefix)_$(flag)_hf_$(νstr)_lk$(lk).jld2")

# ----------------- Hartree-Fock dispersion part ---------------- # 
hf = load(hf_path,"hf");
println("HF energy: ", load(hf_path,"iter_energy")[end])
println("HF convergence: ",load(hf_path,"iter_err")[end])
kvec = reshape(latt.kvec ./ abs(params.g1),lk,lk)
# kvec = reshape(hf.latt.k1,:,1) .+ 1im*reshape(hf.latt.k2,1,:) 
ϵ0 = reshape(hf.ϵk,hf.nt,lk,lk)
# plot_contour_maps(kvec,ϵ0[2,:,:],points=ComplexF64[0+0im],contourlines=[hf.μ],limits=Float64[])
iΓ = (lk%2==0) ? (lk÷2) : ((lk-1)÷2+1)
kcut = real(kvec[:,iΓ])
Ecut = [ϵ0[j,i,iΓ] for j in 1:size(ϵ0,1),i in 1:size(ϵ0,3)];
plot_energy_cuts(kcut,Ecut,lines=[hf.μ])

# ----------------- all the chemical potentials ---------------- # 
νs = collect(-3.9:0.1:3.9)
lk = 20
phi = 0
strain = 2
μs = Float64[]
actual_νs = Float64[]
for ν in νs 
    νstr = round(Int,1000*ν)
    hf_path = joinpath(fpath,"feldman/B0/data/strain$(strain)/phi$(phi)/1_random_hf_$(νstr)_lk$(lk).jld2")
    if ispath(hf_path)
        hf = load(hf_path,"hf");
        push!(actual_νs,round(Int,(hf.ν+4)/8*size(hf.H,1)*size(hf.H,3))/(size(hf.H,1)*size(hf.H,3))*8 - 4)
        push!(μs,hf.μ)
    end
end
fig = figure(figsize=(4,3))
for i in -4:4 
    axvline(i,ls=":",c="gray")
end
axhline(0,ls=":",c="gray")
plot(sort(actual_νs),μs[sortperm(actual_νs)],"b-",ms=2)
xlabel("ν")
ylabel("μ (meV)")
axhline(0,ls="--",c="gray")
# ylim([0,30])
# ylim([0,14])
xlim([-4.2,4.2])
tight_layout()
savefig("cascade_strain_phi$(phi).pdf")
display(fig)
close(fig)