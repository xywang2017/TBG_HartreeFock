using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
# Hartree Fock related 
params = Params(w1=96.056,w0=0.7*96.056,vf=2135.4,dθ=1.05π/180)
w0 = "07"

# ϕs = [1//2;2//5;1//3;2//7;1//4;1//5;1//6;1//8]
ϕs = [1//2;1//3;1//4;1//5;1//6;1//8]
s,t = 3,1

# ------------------------------------------ # 
metadatas = String[]
for ϕ in ϕs 
    flag, seed = "flavor", 1
    if denominator(ϕ)== 4
        flag = "flavor"
    end
    p,q = numerator(ϕ), denominator(ϕ)
    νstr = round(Int,1000*(s+t*p/q))
    metadata = joinpath(fpath,"princeton/data_w$(w0)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    # if !isfile(metadata)
    #     flag = "random"
    #     metadata = joinpath(fpath,"princeton/data_w$(w0)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
    # end
    push!(metadatas,metadata)
end
Δs= plot_spectra_collective(metadatas;savename="spectrum_s$(s)_t$(t).png",titlestr="(s,t)=($(s),$(t))")
