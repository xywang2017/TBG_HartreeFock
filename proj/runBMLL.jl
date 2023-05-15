using PyPlot
using JLD2
fpath = pwd()
include(joinpath(fpath,"libs/bmLL.jl"))

BLAS.set_num_threads(1)

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
    if !isdir(joinpath(fpath,"feldman/B/data_w$(w0str)/_$(p)_$(q)"))
        mkpath(joinpath(fpath,"feldman/B/data_w$(w0str)/_$(p)_$(q)"))
    end
    bm = bmLL()
    nq = (denominator(ϕ)>6) ? 1 : 2
    if q == 3 
        nq = 4 
    elseif q ==2 
        nq = 6
    end
    # if q ==1
    #     nq = 20
    # elseif q ==2 
    #     nq = 20 
    # elseif q==3 
    #     nq = 16 
    # elseif q ==4 
    #     nq = 10 
    # end
    # nq = 16÷q
    println("p= ",p,", q= ",q,", nq= ",nq)
    fname = joinpath(fpath,"feldman/B/data_w$(w0str)/_$(p)_$(q)/_$(p)_$(q)_$(str)_metadata.jld2")
    println(fname)
    # params = Params(w1=96.056,w0=w0*96.056,vf=2135.4,dθ=1.05π/180)
    params = Params(ϵ=0.002,Da=-4100,φ=120.0*π/180,dθ=1.05π/180,w1=110,w0=110*w0,vf=2482)
    initParamsWithStrain(params)
    constructbmLL(bm,params;ϕ= ϕ,nLL=25*q÷p,nq=nq,fname=fname,α=w0, 
        _hBN=false,_strain=true, _σrotation=false, _valley=str,_calculate_overlap=false)
    return bm
end

#
bm = compute_bmLL(ϕ,str,w0,w0str);

#

# jldopen(joinpath(fpath,"feldman/B/data_w07/_1_8/_1_8_K_metadata.jld2")) do file 
    # @time begin 
    #     for m in -3:3, n in -36:36 
    #     Λ = file["$(m)_$(n)"]
    #     F = qr(Λ)
    #     # println(m," ",n," ",norm(Λ))
    #     end
    # end
    # Λ = file["-3_24"]
    # fig = figure(figsize=(5,4))
    # pl=imshow(abs.(Λ),origin="lower")
    # colorbar(pl)
    # axis("equal")
    # display(fig)
    # close(fig)
    # println(norm(Λ))

    # F = svd(Λ)
    # fig = figure()
    # # plot(F.S,"r.")
    # # yscale("log")
    # display(fig)
    # close(fig)
#     energies = file["E"]
#     fig = figure(figsize=(5,4))
#     plot(ones(length(energies)),energies[:],"b_")
#     axis("equal")
#     savefig("test.png",dpi=400)
#     display(fig)
#     close(fig)
# end

# plot spectrum 
function plot_LL_spectrum(ϕs::Vector{Rational{Int}},str::String)
    fig = figure(figsize=(4,3))
    for ϕ in ϕs
        p,q = numerator(ϕ), denominator(ϕ)
        fname = joinpath(fpath,"feldman/B/data_w$(str)/_$(p)_$(q)/_$(p)_$(q)_K_metadata.jld2")
        jldopen(fname) do file 
            energies = file["E"][:]
            tmp_z = real(file["PΣz"])
            chern = zeros(Float64,size(tmp_z,2),size(tmp_z,3),size(tmp_z,4))
            for ik in 1:size(chern,1)
                chern[ik,:,:] = tmp_z[ik,ik,:,:]
            end
            # scatter(ones(length(energies))*ϕ,energies,c=chern[:],cmap="coolwarm",s=2,vmin=-1,vmax=1)
            plot(ones(length(energies))*ϕ,energies,"b.",ms=1)
        end
        # fname = joinpath(fpath,"feldman/B/nonint/data_w$(str)/_$(p)_$(q)/_$(p)_$(q)_Kprime_metadata.jld2")
        # jldopen(fname) do file 
        #     energies = file["E"][:]
        #     # plot(ones(length(energies))*ϕ,energies,"m.",ms=1)
        # end
    end
    lk = 19
    bm_path = joinpath(fpath,"feldman/B0/data/strain2/phi60/bm_lk$(lk).jld2")
    tmp = reshape(load(bm_path,"E"),8,8,lk,lk)
    ϵ0 = zeros(Float64,8,lk,lk)
    for i2 in 1:lk, i1 in 1:lk 
        ϵ0[:,i1,i2] = eigvals(Hermitian(tmp[:,:,i1,i2]))
    end 
    plot(ones(length(ϵ0))*0.0,ϵ0[:],"k.",ms=1)

    ylim([-10,10])
    xlabel(L"ϕ/ϕ_0")
    ylabel("E (meV)")
    tight_layout() 
    savefig("BMLL.pdf")
    display(fig)
    close(fig)
    return nothing
end

plot_LL_spectrum(1 .// collect(1:10) ,"07")
3
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