using PyPlot
using JLD2
fpath = pwd()
include(joinpath(fpath,"libs/bmLL.jl"))

BLAS.set_num_threads(1)

ϕmin = 1//40
str = "K"
w0 = 0.7
w0str = "07"

ϕs = unique(sort([p//q for q in 1:40 for p in 1:q]))
# ϕs = ϕs[ϕs .> ϕmin]
ϕs = ϕs[ϕs .<=0.5]

twist_angle =  138
# calculate spectrum
function compute_bmLL(ϕ::Rational,str::String,w0::Float64,w0str::String)
    fname = "NonInt/Hofstadter/$(twist_angle)_nostrain"
    if !isdir(fname)
        mkpath(fname)
    end
    p = numerator(ϕ)
    q = denominator(ϕ)
    bm = bmLL()
    nq = 40÷denominator(ϕ) 
    # if q ==7 
    #     nq =2  
    # end
    println("p= ",p,", q= ",q,", nq= ",nq)
    # fname = joinpath(fpath,"$(fname)/_$(p)_$(q)_$(str)_metadata.jld2")
    fname = ""
    params = Params(ϵ=0.00,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=110*w0,vf=2482)
    initParamsWithStrain(params)
    constructbmLL(bm,params;ϕ= ϕ,nLL=25*q÷p,nq=nq,fname=fname,α=w0, 
        _hBN=false,_strain=false, _σrotation=false, _valley=str,_calculate_overlap=false)
    return bm.spectrum
end

#
data = Dict()
for ϕ in ϕs
    @time begin
        # println("ϕ: ",ϕ)
        tmp = compute_bmLL(ϕ,str,w0,w0str);
        data["$(ϕ)"] = tmp[:]
    end
end

fname = joinpath(fpath,"NonInt/Hofstadter/$(twist_angle)_nostrain/$(str)_NonIntHofstadter_metadata.jld2")
jldopen(fname, "w") do file
    file["hoftstadter_data"] = data
end


# plot spectrum 
function plot_LL_spectrum()
    fname = joinpath(fpath,"NonInt/Hofstadter/$(twist_angle)_nostrain/K_NonIntHofstadter_metadata.jld2")
    data = load(fname,"hoftstadter_data");
    fig = figure(figsize=(2.5,2.5))
    ϕmin = 1//40
    ϕs = unique(sort([p//q for q in 1:12 for p in 1:q]))
    ϕs = ϕs[ϕs .>= ϕmin]
    ϕs = ϕs[ϕs .<=0.5]
    params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=110*w0,vf=2482)
    initParamsWithStrain(params)
    μB = ZeemanUnit(params)

    for ϕ in ϕs
        p,q = numerator(ϕ), denominator(ϕ)
        energies = reshape(data["$(ϕ)"],2q,:)
        plot(energies[:].-μB*ϕ,ones(length(energies[:]))*ϕ,".",ms=1,markeredgecolor="none",color="tab:red")
        plot(energies[:].+μB*ϕ,ones(length(energies[:]))*ϕ,".",ms=1,markeredgecolor="none",color="tab:blue")
        # plot(ones(length(energies[1:(q-p),:]))*ϕ,energies[1:(q-p),:][:],".",ms=3,markeredgecolor="none",color="r")
        # plot(ones(length(energies[(q-p+1):end,:]))*ϕ,energies[(q-p+1):end,:][:],".",ms=3,markeredgecolor="none",color="gray")
    end
    ylim([0.06,0.51])
    # ylim([-50,50])
    ylabel(L"ϕ/ϕ_0")
    xlabel("E (meV)")
    tight_layout() 
    savefig("test.png",dpi=600,transparent=true)
    display(fig)
    close(fig)
    return nothing
end

plot_LL_spectrum()

# Wannier plot
function plot_wannier(flag=false)
    fname = joinpath(fpath,"NonInt/Hofstadter/$(twist_angle)_nostrain/K_NonIntHofstadter_metadata.jld2")
    data = load(fname,"hoftstadter_data");
    ϕmin = 1//40
    ϕs = unique(sort([p//q for q in 1:40 for p in 1:q]))
    ϕs = ϕs[ϕs .<= 0.5]
    
    params = Params(ϵ=0.00,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=110*w0,vf=2482)
    initParamsWithStrain(params)
    μB = ZeemanUnit(params)

    fig = figure(figsize=(2.5,2.5))
    γ = 0.1   # meV
    E = collect(-50:0.1:50)
    νs = zeros(Float64,length(E),length(ϕs))
    ρνϕ = zeros(Float64,size(νs,1),length(ϕs))
    for iϕ in eachindex(ϕs)
        ϕ = ϕs[iϕ]
        ϵ = data["$(ϕ)"]   
        ϵ = [ϵ .- μB*ϕ; ϵ .+ μB *ϕ]
        for iE in eachindex(E)
            ρνϕ[iE,iϕ] = 1/π * sum(γ./((ϵ .- E[iE]).^2 .+ γ^2))  / (length(ϵ))
            νs[iE,iϕ] = 8/π * sum(atan.((ϵ .- E[iE]) ./ γ)) / (length(ϵ))
        end
    end
    pcolor(νs,repeat(ϕs',size(νs,1)),1 ./ (ρνϕ).^(1/6),cmap="Blues_r" )
    # pcolormesh(νs.*4,ϕs,ρνϕ')
    # colorbar()
    xlim([-4,4])
    xticks(collect(-3:3))
    ylim([0.06,0.51])
    xlabel(L"n/n_s")
    ylabel(L"ϕ/ϕ_0")
    tight_layout()
    display(fig)
    if (flag ==true)
        fname = joinpath(fpath,"Wannier_$(twist_angle)_nostrain.png")
        savefig(fname,dpi=500,transparent=true)
    end
    close(fig)
    return nothing
end

plot_wannier(true)