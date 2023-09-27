using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
# Info and folder name
# ------------------------------------------------------------------------------ # 
twist_angle = 132
foldername = "zeeman/$(twist_angle)_strain"
params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)


w0 = "07"

# ϕs = [1//8;1//6;1//5;1//4;2//7;1//3;2//5;3//7;1//2]
ϕs = sort(unique([p//q for q in 1:12 for p in 1:q]))
ϕs = ϕs[ϕs.<=0.5]
ϕs = ϕs[2:end]
s,t = -3, -1 

# ------------------------------------------------------------------------------ # 

metadatas = String[]
for ϕ in ϕs 
    p,q = numerator(ϕ), denominator(ϕ)
    νstr = round(Int,1000*(s+t*p/q))
    if abs(s+t*p/q) < 4
        test_str = "nu_$(νstr)"
        if isdir("$(foldername)/_$(p)_$(q)")
            files = readdir("$(foldername)/_$(p)_$(q)")
            f = String[]
            for i in files 
                if occursin(test_str,i)
                    push!(f,"$(foldername)/_$(p)_$(q)"*"/"*i)
                end
            end
            println(f)
        end
        # metadata = joinpath(fpath,"$(foldername)/_$(p)_$(q)/1_random_tL_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
        # if !isfile(metadata)
        #     metadata = joinpath(fpath,"$(foldername)/_$(p)_$(q)/1_bm_cascade_tL_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
        # end
        # if isfile(metadata)
        #     E = load(metadata,"iter_energy")[end]
        #     for flag in ["flavor","random","chern","bm","strong","bm_cascade_tL"], seed in 1:10
        #         metadata0 = joinpath(fpath,"$(foldername)/_$(p)_$(q)/$(seed)_$(flag)_init_HF_$(p)_$(q)_nu_$(νstr).jld2")
        #         if isfile(metadata0)
        #             E0 = load(metadata0,"iter_energy")[end]
        #             if E0<=E 
        #                 E, metadata = E0, metadata0 
        #             end
        #         end
        #     end
        #     if load(metadata,"iter_err")[end] > 1e-6
        #         # println("s= ",s," t=",t," p=",p," q=",q," Iter err: ",load(metadata,"iter_err")[end])
        #     end
        #     # println(metadata)
        #     push!(metadatas,metadata)
        # end
    end
end
# Δs= plot_spectra_collective(metadatas;savename="spectrum_s$(s)_t$(t).png",titlestr="(s,t)=($(s),$(t))");