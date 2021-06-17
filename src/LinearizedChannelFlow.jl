include(srcdir("SpectralDiff.jl"))
using LinearAlgebra
using ApproxFun, SparseArrays

function chflow_operator(Nx::Integer,Ny::Integer,Nz::Integer,R::Real)
    N = Nx*(Ny-2)*Nz
    y = cheb_grid(Ny)
    U = (1 .-y.^2)[2:Ny-1]
    Ix = sparse(I,Nx,Nx)
    Iy = sparse(I,Ny-2,Ny-2)
    Iz = sparse(I,Nz,Nz)

    Dx = sparse(Derivative(Fourier(),1)[1:Nx,1:Nx])
    Dxx = sparse(Derivative(Fourier(),2)[1:Nx,1:Nx])
    Dxxxx = sparse(Derivative(Fourier(),4)[1:Nx,1:Nx])
    ∂x = kron(Iz,Iy,Dx)
    ∂xx = kron(Iz,Iy,Dxx)
    ∂xxxx = kron(Iz,Iy,Dxxxx)

    Dy = cheb_diff(Ny)
    Dyy = Dy^2
    s = Diagonal([0; 1 ./(1 .-y[2:Ny-1].^2); 0])
    Dyyyy = (Diagonal(1 .- y.^2)*Dy^4 - 8*Diagonal(y)*Dy^3 - 12*Dy^2)*s
    ∂yy = kron(Iz,Dyy[2:Ny-1,2:Ny-1],Ix)
    ∂yyyy = kron(Iz,Dyyyy[2:Ny-1,2:Ny-1],Ix)

    Dz = sparse(Derivative(Fourier(),2)[1:Nz,1:Nz])
    Dzz = sparse(Derivative(Fourier(),2)[1:Nz,1:Nz])
    Dzzzz = sparse(Derivative(Fourier(),4)[1:Nz,1:Nz])
    ∂z = kron(Dz,Iy,Ix)
    ∂zz = kron(Dzz,Iy,Ix)
    ∂zzzz = kron(Dzzzz,Iy,Ix)

    Δ = ∂xx + ∂yy + ∂zz
    Δ2 = ∂xxxx+∂yyyy+∂zzzz+2*(∂xx*∂yy+∂xx*∂zz+∂yy*∂zz)
    𝒰 = kron(Iz,Diagonal(U),Ix)
    𝒰p = kron(Iz,Diagonal(Dy[2:Ny-1,2:Ny-1]*U),Ix)
    𝒰pp = kron(Iz,Diagonal(Dyy[2:Ny-1,2:Ny-1]*U),Ix)
    𝒪 = spzeros(N,N)

    blocks =[]
    Nb = Nx*(Ny-2)
    for k=1:Nz
        push!(blocks, sparse(inv(Array(Δ[1+(k-1)*Nb:k*Nb,1+(k-1)*Nb:k*Nb]))))
    end
    iΔ = blockdiag(blocks...)

    ℒos = 𝒰*∂x*Δ-𝒰pp*∂x-1/R*Δ2
    ℒsq = -𝒰*∂x+1/R*Δ
    ℒ = [-iΔ*ℒos 𝒪;
         -𝒰p*∂z ℒsq]
    return ℒ
end
