using PyPlot
using JLD2
fpath = joinpath(pwd(),"B0")
include(joinpath(fpath,"libs/HF_mod.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))

# ------------------ Specification ------------------ #
lk = 15
params = Params(ϵ=0.00,Da=0,dθ=1.06π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)
latt = Lattice()
initLattice(latt,params;lk=lk)

bm_path = joinpath(fpath,"data/bm.jld2")

# ------------------ non-interacting part ------------------ #
function compute_bm(latt::Lattice,params::Params;fname::String="placeholder.txt")
    blk = HBM()
    initHBM(blk,latt,params;
            lg=9,_σrotation=false,_calculate_overlap=true,_flag_valley="Both",fname=fname)

    kvec = reshape(latt.kvec ./ abs(params.g1),lk,lk)
    ϵ0 = reshape(blk.hbm,blk.nη,blk.nb,lk,lk)
    # plot_contour_maps(kvec,ϵ0[2,2,:,:],points=[params.Kt/abs(params.g1)])
    iΓ = (lk%2==0) ? (lk÷2) : ((lk-1)÷2+1)
    kcut = real(kvec[:,iΓ])
    Ecut = reshape(ϵ0[:,:,:,iΓ],:,length(kcut))
    plot_energy_cuts(kcut,Ecut)
    return nothing
end

compute_bm(latt,params,fname=bm_path)

# ------------------ Hartree-Fock part ------------------ #
function compute_hf(ν::Float64,latt::Lattice,params::Params;fname::String="placeholder.txt")
    hf = HartreeFock()
    run_HartreeFock(hf,params,latt,bm_path; ν=ν,_Init="random")

    # kvec = reshape(latt.kvec ./ abs(params.g1),lk,lk)
    # ϵ0 = reshape(hf.ϵk,hf.nt,lk,lk)
    # # plot_contour_maps(kvec,ϵ0[1,:,:],points=[params.Kt/abs(params.g1)])
    # iΓ = (lk%2==0) ? (lk÷2) : ((lk-1)÷2+1)
    # kcut = real(kvec[:,iΓ])
    # Ecut = ϵ0[:,:,iΓ]
    # plot_energy_cuts(kcut,Ecut)

    return hf
end

compute_hf(0.0,latt,params)

# ----------------- Hartree-Fock analysis part ---------------- # 
