using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
# Info and folder name
# ------------------------------------------------------------------------------ # 
twist_angle = 120
foldername = "zeeman/$(twist_angle)_strain"
params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=twist_angle*0.01*π/180,w1=110,w0=77,vf=2482)
initParamsWithStrain(params)


w0 = "07"

ϕs = sort(unique([p//q for q in 1:12 for p in 1:q]))
ϕs = ϕs[ϕs.<=0.5]
ϕs = ϕs[2:end]
s,t = -3, -1 

# ------------------------------------------------------------------------------ # 
Φs = Float64[]
Ws = Float64[] 
Ns = Float64[]
μs = Float64[]
for ϕ in ϕs 
    p,q = numerator(ϕ), denominator(ϕ)
    νstr = round(Int,1000*(s+t*p/q))
    metadata = find_lowest_energy_datafile("$(foldername)/_$(p)_$(q)";test_str="nu_$(νstr)",_printinfo=false)
    if !isempty(metadata)
        push!(Ws,load(metadata,"iter_energy")[end])
        push!(Φs,ϕ)
        push!(Ns,(s+t*p/q))
        push!(μs,load(metadata,"hf").μ)
    end
end


# Bohr magneton is μB=eħ/2me [J/T]
# E0 = ħ^2/2me area_moire
aa = 2.46e-10
ħ = 1.054571817e-34
ee = 1.602176634e-19
me = 9.1093837e-31
E0 = ħ^2/(2me*aa^2*params.area) /ee * 1000 # [meV]

fig = figure(figsize=(4,3))
ϕavg = (Φs[1:(end-1)].+Φs[2:end])/2
navg = s .+ t*ϕavg
α = -diff(Ws-μs[end-8]*Ns)./diff(Φs)
plot(ϕavg,α /(2π*E0) ./navg,
            "^-",c="b",ms=3,markeredgecolor="none")
xlabel(L"\rm ϕ/ϕ_0")
ylabel(L"\rm {\bf M}\ (μ_B/e)")
tight_layout()
display(fig)
close(fig)