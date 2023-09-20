using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
# Hartree Fock related 
params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=1.05π/180,w1=110,w0=110*0.7,vf=2482)
initParamsWithStrain(params)
p, q = 1,8
jldopen(joinpath(fpath,"NonInt/124_strain/_$(p)_$(q)_K_metadata.jld2")) do file 
    m,n = 1,-3
    r1 = 4
    Λ = file["$(m)_$(n)"]
    fig = figure(figsize=(5,4))
    pl=imshow(abs.(Λ),origin="lower")
    G = abs(params.g1*m+params.g2*n/q)
    colorbar(pl)
    axis("equal")
    display(fig)
    close(fig)

    tmp = reshape(Λ,2q,q,2q,q);
    lol = tmp[:,4+r1,:,2+r1] ./ (tmp[:,4,:,2]*exp(-1im*2π*r1*n/q))
    # println(sum(real(lol))/(2q)^2)
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