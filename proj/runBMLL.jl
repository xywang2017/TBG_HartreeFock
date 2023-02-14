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
    if !isdir(joinpath(fpath,"feldman/data_w$(w0str)/_$(p)_$(q)"))
        mkpath(joinpath(fpath,"feldman/data_w$(w0str)/_$(p)_$(q)"))
    end
    bm = bmLL()
    nq = (denominator(ϕ)>10) ? 2 : 2
    if denominator(ϕ)<4
        nq = 20 
    end
    if denominator(ϕ)<3
        nq = 40 
    end
    println("p= ",p,", q= ",q,", nq= ",nq)
    fname = joinpath(fpath,"feldman/data_w$(w0str)/_$(p)_$(q)/_$(p)_$(q)_$(str)_metadata.jld2")
    # fname = joinpath(fpath,"feldman/nonint/data_w$(w0str)/_$(p)_$(q)/_$(p)_$(q)_$(str)_metadata.jld2")
    println(fname)
    params = Params(w1=110.0,w0=110.0*w0)
    initParamsWithStrain(params)
    constructbmLL(bm,params;ϕ= ϕ,nLL=25*q÷p,nq=nq,fname=fname,α=w0, 
        _hBN=false,_strain=true, _σrotation=false, _valley=str,_calculate_overlap=true)
    return bm
end

#
bm = compute_bmLL(ϕ,str,w0,w0str);

##

# jldopen(joinpath(fpath,"feldman/data_w07/_1_4/_1_4_K_metadata.jld2")) do file 
#     # for m in -3:3, n in -12:12 
#     #     Λ = file["$(m)_$(n)"]
#     #     if n%4 !=0
#     #         println(norm(tr(Λ)))
#     #     end
#     # end
#     Λ = file["0_0"]
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

# # plot spectrum 
# function plot_LL_spectrum(ϕs::Vector{Rational{Int}},str::String)
#     fig = figure(figsize=(6,4))
#     for ϕ in ϕs
#         p,q = numerator(ϕ), denominator(ϕ)
#         fname = joinpath(fpath,"feldman/nonint/data_w$(str)/_$(p)_$(q)/_$(p)_$(q)_K_metadata.jld2")
#         jldopen(fname) do file 
#             energies = file["E"][:]
#             plot(ones(length(energies))*ϕ,energies,"g.",ms=1)
#         end
#         # fname = joinpath(fpath,"feldman/nonint/data_w$(str)/_$(p)_$(q)/_$(p)_$(q)_Kprime_metadata.jld2")
#         # jldopen(fname) do file 
#         #     energies = file["E"][:]
#         #     # plot(ones(length(energies))*ϕ,energies,"m.",ms=1)
#         # end
#     end
#     xlabel(L"ϕ/ϕ_0")
#     ylabel("E (meV)")
#     tight_layout() 
#     savefig("feldman/figures/BM_spectrum.pdf")
#     display(fig)
#     close(fig)
#     return nothing
# end

# plot_LL_spectrum(1 .// collect(1:17),"06")


# ## weak coupling (0,2) gaps 
# function find_gaps(ϕs::Vector{Rational{Int}},str::String)
#     s=0
#     fig = figure(figsize=(4,3))
#     ts = [0,1,2]
#     for t in ts 
#         gaps = Float64[]
#         for ϕ in ϕs
#             p,q = numerator(ϕ), denominator(ϕ)
#             fnameK = joinpath(fpath,"yacoby/nonint/data_w$(str)/_$(p)_$(q)/_$(p)_$(q)_K_metadata.jld2")
#             fnameKprime = joinpath(fpath,"yacoby/nonint/data_w$(str)/_$(p)_$(q)/_$(p)_$(q)_Kprime_metadata.jld2")
#             energies = Float64[]
#             jldopen(fnameK) do file 
#                 energies = [energies ;file["E"][:]]
#             end
#             jldopen(fnameKprime) do file 
#                 energies = [energies ;file["E"][:]]
#             end
#             energies = sort(energies)
#             Nν = round(Int, ( 4 + (s+t*ϕ) ) / 8 * length(energies) ) 
#             push!(gaps,energies[Nν+1]-energies[Nν])
#         end
#         plot(ϕs,gaps,"-",label="($s,$t)")
#     end
#     legend()
#     xlabel(L"ϕ/ϕ_0")
#     ylabel("Δ (meV)")
#     # ylim([0,2.5])
#     tight_layout() 
#     savefig("ratio_0665_w0_08.pdf")
#     display(fig)
#     close(fig)
#     return nothing
# end

# find_gaps(1 .// collect(2:20),"08")