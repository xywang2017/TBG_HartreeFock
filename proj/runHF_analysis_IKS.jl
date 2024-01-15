using PyPlot,JLD2
using Printf 
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))

#
# Info and folder name
# ------------------------------------------------------------------------------ # 
twist_angles = [105; collect(106:2:138)] 
twist_angle = 120
# for twist_angle in twist_angles
# dir = "/media/xiaoyuw@ad.magnet.fsu.edu/Data/Code/TBG_HartreeFock/"
dir = "/Volumes/Data/Code/TBG_HartreeFock/"
dir = ""
foldername = dir*"MinHao/$(twist_angle)_strain"
params = Params(ϵ=0.001,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)

# ----------------------------------Hartree Fock spectrum-------------------------------------------- # 
s,t = -2,-2
p,q = 1,3
νF = (s)+(t)*p/q
νstr = round(Int,1000*νF)
metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="random_init_HF_$(p)_$(q)_nu_$(νstr)",_printinfo=true)

plot_density_matrix_bm(metadata,ik=1)
test_tL2_breaking(metadata)
plot_density_matrix_global_order_parameters(metadata)

# ----------------------------------IKS Analysis-------------------------------------------- # 

P = reshape(load(metadata,"hf").P,2q,2,2,2q,2,2,:);

p_subblock = P[:,2,1,:,2,1,:];

_pratio = p_subblock[:,:,4]./ p_subblock[:,:,6]
