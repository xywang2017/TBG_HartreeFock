# using PyPlot
using JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))

str = ARGS[1]
w0 = parse(Float64,ARGS[2])*0.1
w0str = ARGS[2]
p = parse(Int,ARGS[3])
q = parse(Int,ARGS[4])
ϕ = p//q

# calculate spectrum
function compute_bmLL(ϕ::Rational,str::String,w0::Float64,w0str::String)
    p = numerator(ϕ)
    q = denominator(ϕ)
    if !isdir(joinpath(fpath,"data_w$(w0str)/_$(p)_$(q)"))
        mkdir(joinpath(fpath,"data_w$(w0str)/_$(p)_$(q)"))
    end
    bm = bmLL()
    nq = (denominator(ϕ)>6) ? 1 : 1
    println("p= ",p,", q= ",q,", nq= ",nq)
    fname = joinpath(fpath,"data_w$(w0str)/_$(p)_$(q)/_$(p)_$(q)_$(str)_metadata.jld2")
    println(fname)
    constructbmLL(bm;ϕ= ϕ,nLL=25*q÷p,nq=nq,fname=fname,α=w0, 
        _hBN=false,_strain=false, _σrotation=false, _valley=str,_calculate_overlap=true)
    return bm
end

#
bm = compute_bmLL(ϕ,str,w0,w0str);

# ##

# jldopen(joinpath(fpath,"data_w02/_1_7/_1_7_K_metadata.jld2")) do file 
#     # for m in -3:3, n in -12:12 
#     #     Λ = file["$(m)_$(n)"]
#     #     if n%4 !=0
#     #         println(norm(tr(Λ)))
#     #     end
#     # end
#     Λ = file["3_21"]
#     fig = figure(figsize=(5,4))
#     pl=imshow(abs.(Λ),origin="lower")
#     colorbar(pl)
#     axis("equal")
#     display(fig)
#     close(fig)
#     println(norm(Λ))

#     energies = file["E"]
#     fig = figure(figsize=(5,4))
#     plot(ones(length(energies)),energies[:],"b_")
#     colorbar(pl)
#     axis("equal")
#     display(fig)
#     close(fig)
# end
