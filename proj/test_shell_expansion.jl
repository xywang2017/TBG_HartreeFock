using PyPlot,JLD2
fpath = pwd()
include(joinpath(fpath,"libs/MagneticFieldHF.jl"))
include(joinpath(fpath,"libs/plot_helpers.jl"))
#
# Hartree Fock related 
params = Params(ϵ=0.002,Da=-4100,φ=0.0*π/180,dθ=1.20π/180,w1=110,w0=110*0.7,vf=2482)
initParamsWithStrain(params)
p, q = 2,7
jldopen(joinpath(fpath,"NonInt/120_strain/_$(p)_$(q)_Kprime_metadata.jld2")) do file 
    m,n = 3,-3q
    r1 = 3
    Λ = file["$(m)_$(n)"]
    # Λ = file["PΣz"]
    fig = figure(figsize=(5,4))
    pl=imshow(abs.(Λ),origin="lower")
    G = abs(params.g1*m+params.g2*n/q)
    colorbar(pl)
    axis("equal")
    display(fig)
    close(fig)
end

for m in -3:3, n in -15:15
    Λ1 = load(joinpath(fpath,"NonInt/120_strain/_$(p)_$(q)_Kprime_test_metadata_v1.jld2"),"$(m)_$(n)");
    Λ2 = load(joinpath(fpath,"NonInt/120_strain/_$(p)_$(q)_Kprime_metadata_v1.jld2"),"$(m)_$(n)");
    
    Λ1 = reshape(Λ1,2q,q,4,2q,q,4)[:,:,2,:,:,1]
    Λ2 = reshape(Λ2,2q,q,4,2q,q,4)[:,:,2,:,:,1]

    for j in 1:2q, i in 1:2q
        θ = Λ1[i,3,j,2] / Λ2[i,3,j,2]
        tmp = Λ1[i,:,j,:] -  Λ2[i,:,j,:]*θ
        if norm(tmp)> 1e-6
            println(m," ",n," ",i," ",j," error with C2P: ",norm(tmp))
            # println(m," ",n)
        end
    end
    # display("text/plain",Λ1[1,:,5,:]./ Λ2[1,:,5,:])
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