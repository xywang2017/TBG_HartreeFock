using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
# Hartree Fock related 
params = Params(w1=96.056,w0=0.7*96.056,vf=2135.4,dθ=1.05π/180)

##
flag = "random"
seed = 4
w0s = ["07"]
w0snum = [0.7]
p,q = 1,8
νF = 3 +(-1)*p/q
νstr = round(Int,1000*νF)
hf = 0
for w0 in w0s
    metadata = joinpath(fpath,"princeton/data_w$(w0)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    println("HF energy: ",load(metadata,"iter_energy")[end])
    println("Convergence: ",load(metadata,"iter_err")[end])
    # hf = load(metadata,"hf")
    # iter_err = load(metadata,"iter_err")
    # iter_energy = load(metadata,"iter_energy")
    # iter_oda = load(metadata,"iter_oda")
    # hf.Λs =  Array{ComplexF64,4}(undef,0,0,0,0)
    # hf.Λ =  Array{ComplexF64,4}(undef,0,0,0,0)
    # save(metadata,"hf",hf,
    #                 "iter_err",iter_err,"iter_energy",iter_energy,"iter_oda",iter_oda)
    # plot_hf_iterations(metadata)
    plot_spectra(metadata;savename="test.pdf")
    # hf = load(metadata,"hf")
end 

## BM basis 
seed = 1
p, q = 1,5
flag = "flavor"
νF = 3 + (0)*p/q
νstr = round(Int,1000*νF)
metadata = joinpath(fpath,"princeton/data_w07/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
# plot_density_matrix_bm_valley_spin(metadata)
# plot_density_matrix_bm(metadata)
plot_density_matrix_sublattice(metadata)
plot_density_matrix_sublattice_full(metadata)


## strong coupling basis at reference point defined by metadata0

seed = 1
flag  = "rand"
νF0 = 0+ (0)*p/q
νstr0 = round(Int,1000*νF0)
metadata0 = joinpath(fpath,"princeton/data_w07/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr0).jld2")
# plot_density_matrix_strong_coupling(metadata,metadata0)
# plot_density_matrix_strong_coupling_valley_spin(metadata,metadata0)
# 
# plot_order_parameters(metadata)


##
x = [0,1,2,3,4,5,6,7]
y = [13.23;5.05;7.36;8.26;5.92;4.72;16.44;1.92];
y8 = [12;2.57;5.15;8.05;3.84;3.51;12.57;3.09]
fig = figure(figsize=(4,3))
plot(x,y,"b^",label="1/4")
plot(x,y8,"r<",label="1/8")
xticks(x,["(0,0)","(0,1)","(0,2)","(1,-1)","(1,0)","(2,0)","(2,1)","(2,2)"])
ylim([0,19])
yticks(0:4:20)
legend()
xlabel("(s,t)")
ylabel("Δ (meV)")
tight_layout()
savefig("gaps_strongcoupling.pdf")
display(fig)
close(fig)


## 

# ------------------ plot_spectra_collectively_at_different_flux --------------- # 
ϕs = [1//4;1//8]
w0 = "07"
metadatas = String[]
for ϕ in ϕs 
    flag, seed = "flavor", 1
    p,q = numerator(ϕ), denominator(ϕ)
    s,t = 1,0
    νstr = round(Int,1000*(s+t*p/q))
    metadata = joinpath(fpath,"princeton/data_w$(w0)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    push!(metadatas,metadata)
end
Δs= plot_spectra_collective(metadatas;savename="tmp.pdf",titlestr="(s,t)=(2,0)")