mutable struct HybridWannier 
    #|γ;k⟩ = |m,k⟩* Wγm_{m,γ};   Wγk is 2x2xlkxlk matrix 
    #note that Wγm = is properly normalized
    Wγm::Array{ComplexF64,4}    
    Wη::Array{ComplexF64,3}  # Wilson loop along the per k2 value is a 2x2 matrix, 2x2xlk
    αη::Array{ComplexF64,3}  # coefficients
    ϵη::Array{ComplexF64,2}  # the eigenvalue of the Wilson loop associated with the eigenvector [1,i]
    Λη::Array{ComplexF64,4}  # partial Wilson loop
    WγLS::Array{ComplexF64,4} # Wγk in layer-sublattice basis
    Hγk::Array{ComplexF64,4} # the kinetic energy in the chiral basis, in momentum space
    HybridWannier() = new()     
end


function initHybridWannier(bm::HBM,iη::Int)
    """
    Construct hybrid Wannier states from Bloch eigenstates C2T + 1
    for valley iη
    """
    lk = Int(sqrt(length(bm.latt.kvec)))
    l1 = lk
    l2 = lk
    lg = bm.lg
    nlocal = 4

    A = HybridWannier()
    A.Wη = zeros(ComplexF64,2,2,l2)  # Wilson loop
    A.αη = zeros(ComplexF64,2,2,l2)
    A.ϵη = zeros(ComplexF64,2,l2)  # the eigenvalue of the Wilson loop associated with the eigenvector [1,i]
    A.Λη = zeros(ComplexF64,2,2,l1,l2) # partial Wilson loop
    A.Wγm = zeros(ComplexF64,2,2,l1,l2)
    A.WγLS = zeros(ComplexF64,4*lg^2,2,l1,l2) # hybrid wannier states in layer-subbm.lattice basis
    A.Hγk = zeros(ComplexF64,2,2,l1,l2)
    Uk = reshape(view(bm.Uk,:,:,iη,:),:,2,l1,l2)

    for ik in 1:l2 
        A.αη[:,:,ik] = ComplexF64[1 1; 1im -1im] / sqrt(2)
    end
    
    # smooth gauge along g2
    for ik in 1:(l2-1)
        for n in 1:2
            u1 = view(Uk,:,n,1,ik)
            u2 = view(Uk,:,n,1,ik+1)
            if real((u1'*u2)[1,1])<0
                u2 .*= -1
            end
        end
    end

    for ik in 1:(l2-1)
        u1 = view(Uk,:,2,1,ik)
        u2 = view(Uk,:,2,1,ik+1)
        if real((u1'*u2)[1,1])<0
            u2 .*= -1
        end
    end
    
    # cut along g2
    for ik in 1:l2
        A.Λη[:,:,1,ik] .= Array{ComplexF64}(I,2,2)
        for jk in 1:l1 # along g1 direction
            u1 = view(Uk,:,:,jk,ik)
            if (jk<l1)
                u2 = view(Uk,:,:,jk+1,ik)
                F = svd(u1'*u2)
                A.Λη[:,:,jk+1,ik] .= A.Λη[:,:,jk,ik]*(F.U*F.Vt)
            else 
                #(N-1,1)
                u2 = reshape(view(Uk,:,:,1,ik),nlocal,lg,lg,2)
                u2 = reshape(circshift(u2,(0,-1,0,0)), nlocal*lg^2,2)
                F = svd( u1'*u2 )
                A.Wη[:,:,ik] .= A.Λη[:,:,jk,ik] * (F.U*F.Vt)
            end
        end
        F = eigen(A.Wη[:,:,ik])
        vals = A.αη[:,:,ik]' * A.Wη[:,:,ik] * A.αη[:,:,ik]
        if (abs(vals[1,2]+abs(vals[2,1]))>1e-5)
            print("(1,i) is not the eigenvector\n")
            A.αη[:,1,ik] = normalize(F.vectors[:,1])
            A.αη[:,2,ik] = normalize(F.vectors[:,2])
            A.ϵη[1,ik] = F.values[1]
            A.ϵη[2,ik] = F.values[2]
        else
            A.ϵη[1,ik] = vals[1,1]
            A.ϵη[2,ik] = vals[2,2]
        end
    end
    
    # construct wannier states
    for ik in 1:l2
        for j in 1:l1
            A.Wγm[:,:,j,ik] = inv(A.Λη[:,:,j,ik]) * A.αη[:,:,ik] * Diagonal(A.ϵη[:,ik].^((j-1)/l1)) 
            u = view(Uk,:,:,j,ik)
            A.WγLS[:,:,j,ik] = u*A.Wγm[:,:,j,ik]
        end 
        pA = norm(A.WγLS[1:2:(4*lg^2),1,1,ik])
        pB = norm(A.WγLS[2:2:(4*lg^2),1,1,ik])
        # print("\n")
        # print(pA," ",pB,"\n")
        if pA < pB
            A.ϵη[:,ik] = A.ϵη[2:-1:1,ik]
            A.αη[:,:,ik] = A.αη[:,2:-1:1,ik]
            A.Wγm[:,:,:,ik] = A.Wγm[:,2:-1:1,:,ik]
            A.WγLS[:,:,:,ik] = A.WγLS[:,2:-1:1,:,ik]
        end
    end
    
    for i2 in 1:l2
        for i1 in 1:l1
            A.Hγk[:,:,i1,i2] =  A.Wγm[:,:,i1,i2]' * Diagonal(bm.spectrum[:,iη,(i2-1)*l1+i1]) * A.Wγm[:,:,i1,i2]
        end
    end
    return A
end

mutable struct HybridWannierRealSpace 
    ndim::Int  # number of n values to consider
    nvec::Vector{Int}  # set of n values to perform fourier tranform 
    area::Array{ComplexF64,3} # field of view in absolute units
    rdim::Int # total size of the field of view
    area_extent::Int # field of view boundaries, typically 2bm.lattice vector on either side
    Wγn::Array{ComplexF64,6} # 4x2xlk (g2) x ndim 
    HybridWannierRealSpace() = new()
end

@inline function FourierTransform(z::ComplexF64,A::Array{ComplexF64,5},bm::HBM) 
    expikr = reshape(exp.(1im * real.(bm.latt.kvec*z')),1,1,bm.latt.lk,bm.latt.lk)
    expigr = reshape(exp.(1im*real.(bm.gvec*z') ),1,bm.lg^2,1,1,1)
    B = reshape(sum(reshape(sum(expigr .* A,dims=2),4,2,bm.latt.lk,bm.latt.lk) .* expikr, dims=3),4,2,bm.latt.lk)
    return B::Array{ComplexF64,3}
end

function initHybridWannierRealSpace(bm::HBM;area_extent::Int=2,nvec::Vector{Int}=[0],rdim::Int=21)
    lk = bm.latt.lk
    lg = bm.lg 
    params = bm.params

    HW_R = HybridWannierRealSpace()
    HW_R.nvec = nvec
    HW_R.ndim = length(HW_R.nvec)
    HW_R.area_extent = area_extent
    HW_R.rdim = rdim # r-grid along both a1 and a2 directions
    HW_R.area = zeros(ComplexF64,HW_R.rdim,HW_R.rdim,HW_R.ndim)
    HW_R.Wγn = zeros(ComplexF64,4,HW_R.rdim,HW_R.rdim,2,lk,HW_R.ndim)
    for n in 1:size(HW_R.area,3)
        # extent1 = range((-HW_R.area_extent+HW_R.nvec[n])*params.a1,(HW_R.area_extent+HW_R.nvec[n])*params.a1,length=HW_R.rdim)
        extent1 = range((-HW_R.area_extent+HW_R.nvec[n])*real(params.a1),(HW_R.area_extent+HW_R.nvec[n])*real(params.a1),length=HW_R.rdim)
        extent2 = range(-1.5params.a2,1.5params.a2,length=HW_R.rdim)  # periodicity along g2 means only need to keep 1 unit cell
        HW_R.area[:,:,n] = reshape(extent1,HW_R.rdim,1) .+ reshape(extent2,1,HW_R.rdim)
    end
    
    Wγn_Components = zeros(ComplexF64,4,lg^2,2,lk,lk)
    Uk = reshape(bm.Uk,:,2,2,bm.latt.lk,bm.latt.lk)
    for n in 1:size(HW_R.area,3)
        for ik in 1:lk 
            for ik2 in 1:lk
                Wγn_Components[:,:,:,ik,ik2] =  reshape( view(Uk,:,:,1,ik,ik2) * 
                                                         exp(-1im * 2π * (ik-1)* HW_R.nvec[n] /sqrt(lk) ), (4,lg^2,2) )
            end
        end
        for i2 in 1:size(HW_R.area,2)
            for i1 in 1:size(HW_R.area,1)
                z0 = HW_R.area[i1,i2,n]
                HW_R.Wγn[:,i1,i2,:,:,n] = FourierTransform(z0,Wγn_Components,bm) 
            end
        end
    end

    return HW_R
end

function plot_HybridWannier_RealSpace(HW_R::HybridWannierRealSpace,params::Params)
    n = 1
    wpmR00_p = abs2.(view(HW_R.Wγn,:,:,:,1,1,n))
    zmesh = HW_R.area[:,:,n]
    xgrid = real.(zmesh)./abs(params.a1)
    ygrid = imag.(zmesh)./abs(params.a1)

    points_r  = -2:2 
    latt_points = ( reshape(points_r,:,1) * params.a1 .+ reshape(points_r,1,:) * params.a2) ./ abs(params.a1)

    fig, ax = subplots(2,2,figsize=(6,6))
    pl=ax[1,1].contourf(xgrid,ygrid,wpmR00_p[1,:,:])
    ax[1,1].scatter(real(latt_points),imag(latt_points),c="gray",s=3)
    colorbar(pl,ax=ax[1,1])
    ax[1,1].axis("equal")
    ax[1,1].set_title("Top A")
    pl=ax[1,2].contourf(xgrid,ygrid,wpmR00_p[2,:,:])
    ax[1,2].scatter(real(latt_points),imag(latt_points),c="gray",s=3)
    colorbar(pl,ax=ax[1,2])
    ax[1,2].axis("equal")
    ax[1,2].set_title("Top B")
    pl=ax[2,1].contourf(xgrid,ygrid,wpmR00_p[3,:,:])
    ax[2,1].scatter(real(latt_points),imag(latt_points),c="gray",s=3)
    colorbar(pl,ax=ax[2,1])
    ax[2,1].axis("equal")
    ax[2,1].set_title("Bottom A")
    pl=ax[2,2].contourf(xgrid,ygrid,wpmR00_p[4,:,:])
    ax[2,2].scatter(real(latt_points),imag(latt_points),c="gray",s=3)
    colorbar(pl,ax=ax[2,2])
    ax[2,2].axis("equal")
    ax[2,2].set_title("Bottom B")
    tight_layout()
    display(fig)
    # savefig("hybrid_wannier_example_v1.pdf")
    close(fig)
end

function plot_HybridWannier_RealSpaceCombined(HW_R::HybridWannierRealSpace,params::Params;savename="hybrid_wannier_example_v1.pdf",nk::Int=1)
    n = 1
    zmesh = HW_R.area[:,:,n]
    xgrid = real.(zmesh)./abs(params.a1)
    ygrid = imag.(zmesh)./abs(params.a1)
    max_val = maximum(sum(abs2.(view(HW_R.Wγn,:,:,:,:,:,n)),dims=1))
    wpmR00_p = reshape(sum(abs2.(view(HW_R.Wγn,:,:,:,1,nk,n)),dims=1),size(zmesh,1),size(zmesh,2)) 
    wpmR00_m = reshape(sum(abs2.(view(HW_R.Wγn,:,:,:,2,nk,n)),dims=1),size(zmesh,1),size(zmesh,2)) 
    wpmR00_p ./= max_val #maximum(wpmR00_p)
    wpmR00_m ./= max_val
    # wpmR00_p = reshape(abs2.(view(HW_R.Wγn,1,:,:,1,1,n)),size(zmesh,1),size(zmesh,2))

    points_r  = -2:2 
    latt_points = ( reshape(points_r,:,1) * params.a1 .+ reshape(points_r,1,:) * params.a2) ./ abs(params.a1)

    fig,ax = subplots(2,1,figsize=(6,4))
    # contourf(xgrid,ygrid,wpmR00_p,cmap="Blues")
    pl=ax[1,1].pcolormesh(xgrid,ygrid,wpmR00_p,cmap="Blues",shading="gouraud")
    ax[1,1].scatter(real(latt_points),imag(latt_points),c="k",s=3)
    # colorbar(pl,fraction=0.03,location="left")
    # title(L"$|w_{0,\frac{%$(nk)}{64}}(\mathbf{r})|^2$")
    ax[1,1].axis("equal")
    ax[1,1].set_xlim([-0.8,0.6])
    ax[1,1].set_ylim([-1.4,1.4])
    ax[1,1].axis("off")

    pl=ax[2,1].pcolormesh(xgrid,ygrid,wpmR00_m,cmap="Reds",shading="gouraud")
    ax[2,1].scatter(real(latt_points),imag(latt_points),c="k",s=3)
    # colorbar(pl,fraction=0.03,location="left")
    # title(L"$|w_{0,\frac{%$(nk)}{64}}(\mathbf{r})|^2$")
    ax[2,1].axis("equal")
    ax[2,1].set_xlim([-0.8,0.6])
    ax[2,1].set_ylim([-1.4,1.4])
    ax[2,1].axis("off")

    tight_layout()
    display(fig)
    # savefig(savename,dpi=600,transparent=false)
    close(fig)
end