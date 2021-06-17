using DrWatson
@quickactivate "Data-driven resolvent"

include(srcdir("SpectralDiff.jl"))
using ApproxFun, SparseArrays
using KrylovKit
using DifferentialEquations

##
Sf = Fourier()
Sc = Chebyshev()
Sf⊗Sc⊗Sf # ⊗ is \otimes

##

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

function chflow_ode!(dq̂,q̂,ℒ,t)
    mul!(dq̂,ℒ,q̂)
end

function chflow_step(q̂k,ℒ,Δt,tol)::SparseVector{Complex{Float64},Int64}
    prob = ODEProblem{true}(chflow_ode!,q̂k,(0,Δt),ℒ)
    q̂kp1 = similar(q̂k)
    q̂kp1 = solve(prob, RK4(), dt=0.02, save_start=false, saveat=Δt, reltol=tol, abstol=tol).u[1]
    return q̂kp1
end

Nx = 8
Ny = 33
Nz = 8
R = 2000
ℒ = chflow_operator(Nx,Ny,Nz,R)

N = 2*Nx*(Ny-2)*Nz
q̂0 = sparse(randn(Complex{Float64},N))
Δt = 2.0
tol = 1e-6
A(q̂k) = chflow_step(q̂k,ℒ,Δt,tol)
@time A(q̂0)
KrylovDefaults.tol
KrylovDefaults.maxiter
KrylovDefaults.krylovdim
@time λ, V, info = eigsolve(ℒ, q̂0, 20, :LR, tol=1e-6, krylovdim=1024)
# λ2,V2 = eigen(Array(ℒ))
# scatter(log.(ρ)/Δt, aspectratio=1)
scatter(λ, aspectratio=1)
##
function pois_sys(N::Integer, α::Number, β::Number, R::Number)
    k = sqrt(α^2+β^2)
    D = cheb_diff(N)
    D2 = D^2
    D̃ = D[2:N-1,2:N-1]
    D̃2 = D2[2:N-1,2:N-1]
    M = k^2*I-D̃2
    invM = M\I

    A = pois(N,α,β,R)

    B = [im*α*invM*D̃ invM*k^2 im*β*invM*D̃;
         im*β*I zero(D̃) -im*α*I]

    C = [im*α/k^2*D̃ -im*β/k^2*I;
         I zero(D̃);
         im*β/k^2*D̃ im*α/k^2*I]

    sys = ss(A,B,C,0)
    return sys
end


function pois_weights(N::Integer, k::Number)
    D = cheb_diff(N)
    W = Array(Diagonal(cheb_weights(N)))
    D̃ = D[2:N-1,2:N-1]
    W̃ = W[2:N-1,2:N-1]
    Q = [(k^2*W̃ + D̃'*W̃*D̃) zero(D̃);
          zero(D̃) W̃]/(2*k^2)
    Q = (Q+Q')/2
    return Q
end

function pois_weights(N::Integer)
    W = Array(Diagonal(cheb_weights(N)))
    W̃ = W[2:N-1,2:N-1]
    Q = 0.5*cat(W̃,W̃,W̃,dims=(1,2))
    return Q
end
