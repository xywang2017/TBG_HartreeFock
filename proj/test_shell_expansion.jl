using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
# Hartree Fock related 
params = Params(w1=110,w0=77,vf=2482,dθ=1.05π/180)
initParamsWithStrain(params)
p, q = 1,12
jldopen(joinpath(fpath,"feldman/B/data_w07/_$(p)_$(q)/_$(p)_$(q)_Kprime_metadata.jld2")) do file 
    m,n = 3,-3q
    Λ = file["$(m)_$(n)"]
    fig = figure(figsize=(5,4))
    pl=imshow(abs.(Λ),origin="lower")
    G = abs(params.g1*m+params.g2*n/q)
    println(G)
    colorbar(pl)
    axis("equal")
    display(fig)
    close(fig)
end

q = 3
fig = figure(figsize=(4,4)) 
G0 = abs(3*params.g1+params.g2*3)*1.000001
for m in -3:3, n in -3q:3q 
    G = (params.g1*m+params.g2*(n÷ q) ) 
    GG = (params.g1*m+params.g2*(n/q) ) 
    plot([real(GG)],[imag(GG)],"b.")
    if abs(GG) <G0*cos(pi/6)/abs(cos(mod(angle(GG),pi/3)-pi/6))
        plot([real(GG)],[imag(GG)],"r.")
    end
end
axis("equal")
display(fig)
close(fig)