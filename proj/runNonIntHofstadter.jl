using PyPlot
using JLD2
fpath = pwd()
include(joinpath(fpath,"libs/bmLL.jl"))

BLAS.set_num_threads(1)

str = "K"
w0 = 0.7
w0str = "07"
p = parse(Int,ARGS[1])
q = parse(Int,ARGS[2])
ϕ = p//q

# calculate spectrum
function compute_bmLL(ϕ::Rational,str::String,w0::Float64,w0str::String)
    fname = "105_nostrain"
    p = numerator(ϕ)
    q = denominator(ϕ)
    bm = bmLL()
    nq = 16÷denominator(ϕ) 
    if q ==1
        nq = 20
    elseif q ==2 
        nq = 20 
    elseif q==3 
        nq = 16 
    elseif q ==4 
        nq = 10 
    end
    nq = 16÷q
    println("p= ",p,", q= ",q,", nq= ",nq)
    fname = joinpath(fpath,"$(fname)/B/_$(p)_$(q)_$(str)_metadata.jld2")
    println(fname)
    params = Params(ϵ=0.00,Da=-4100,φ=0.0*π/180,dθ=1.05π/180,w1=110,w0=110*w0,vf=2482)
    initParamsWithStrain(params)
    constructbmLL(bm,params;ϕ= ϕ,nLL=40*q÷p,nq=nq,fname=fname,α=w0, 
        _hBN=false,_strain=false, _σrotation=false, _valley=str,_calculate_overlap=false)
    return bm
end

#
bm = compute_bmLL(ϕ,str,w0,w0str);

#

# plot spectrum 
function plot_LL_spectrum(ϕs::Vector{Rational{Int}},str::String)
    foldername = "105_nostrain"
    fig = figure(figsize=(4,3))
    for ϕ in ϕs
        p,q = numerator(ϕ), denominator(ϕ)
        fname = joinpath(fpath,"$(foldername)/B/_$(p)_$(q)_K_metadata.jld2")
        jldopen(fname) do file 
            energies = file["E"][:]
            tmp_z = real(file["PΣz"])
            chern = zeros(Float64,size(tmp_z,2),size(tmp_z,3),size(tmp_z,4))
            for ik in 1:size(chern,1)
                chern[ik,:,:] = tmp_z[ik,ik,:,:]
            end
            scatter(ones(length(energies))*ϕ,energies,c=chern[:],cmap="coolwarm",s=3,vmin=-1,vmax=1)
        end
    end
    xlabel(L"ϕ/ϕ_0")
    ylabel("E (meV)")
    tight_layout() 
    savefig("BMLL.pdf")
    display(fig)
    close(fig)
    return nothing
end

plot_LL_spectrum(collect(1:16) .//16 ,"07")
