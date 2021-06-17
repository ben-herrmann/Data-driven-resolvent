using DrWatson
@quickactivate "Data-driven resolvent"
using LinearAlgebra, ApproxFun, DifferentialEquations
include(srcdir("NonModalStabilityTools.jl"))
include(srcdir("SpectralDiff.jl"))

# Local Framework

function pois_operator(N::Integer, p)
    α, β, R = p
    k = sqrt(α^2+β^2)
    y = cheb_grid(N)
    U = 1 .-y.^2 #base flow
    D = cheb_diff(N)
    D2 = D^2
    S = Diagonal([0; 1 ./(1 .-y[2:N-1].^2); 0])
    D4 = (Diagonal(1 .- y.^2)*D^4 - 8*Diagonal(y)*D^3 - 12*D^2)*S
    W = Diagonal(cheb_weights(N))

    Ũ = U[2:N-1]
    D̃ = D[2:N-1,2:N-1]
    D̃2 = D2[2:N-1,2:N-1]
    D̃4 = D4[2:N-1,2:N-1]
    W̃ = W[2:N-1,2:N-1]
    Q = [(k^2*W̃ + D̃'*W̃*D̃) zero(D̃);
          zero(D̃) W̃]/(2*k^2)

    M = k^2*I-D̃2
    M2 = k^4*I-2*k^2*D̃2+D̃4
    invM = M\I

    Aos = invM*(-im*α*Diagonal(Ũ)*M-im*α*Diagonal(D̃2*Ũ)-(1/R)*M2)
    Asq = -im*α*Diagonal(Ũ)-(1/R)*M
    Ac = Array(-im*β*Diagonal(D̃*Ũ))

    A = [Aos zero(D̃);
         Ac Asq]
    return A, Q
end

function pois_liftup_setup()
    N = 65
    α = 0.25
    β = 2.0
    R = 2000.0
    p = [α, β, R]
    y = cheb_grid(N)
    A,Q = pois_operator(N,p)
    return A,Q,y,p
end

function pois_TS_setup()
    N = 65
    α = 1.0
    β = 0.0
    R = 2000.0
    p = [α, β, R]
    y = cheb_grid(N)
    A,Q = pois_operator(N,p)
    return A,Q,y,p
end

function q2u(q,α,β)
    k = sqrt(α^2+β^2)
    N = div(length(q),2)+2
    D = cheb_diff(N)
    D2 = D^2
    D̃ = D[2:N-1,2:N-1]
    D̃2 = D2[2:N-1,2:N-1]
    M = k^2*I-D̃2
    invM = M\I
    C = [im*α/k^2*D̃ -im*β/k^2*I;
         I zero(D̃);
         im*β/k^2*D̃ im*α/k^2*I]
    u = C*q
    return u
end

function u2q(u,α,β)
    k = sqrt(α^2+β^2)
    N = div(length(u),3)+2
    D = cheb_diff(N)
    D2 = D^2
    D̃ = D[2:N-1,2:N-1]
    D̃2 = D2[2:N-1,2:N-1]
    M = k^2*I-D̃2
    invM = M\I
    B = [im*α*invM*D̃ invM*k^2 im*β*invM*D̃;
         im*β*I zero(D̃) -im*α*I]
    q = B*u
    return q
end

function q̂2u(q̂,α,β,Nx,Nz)
    k = sqrt(α^2+β^2)
    Ny = div(length(q̂),2)+2
    x = fourier_grid(Nx)
    y = cheb_grid(Ny)
    z = fourier_grid(Nz)
    û = q2u(q̂,α,β)
    û = [û[1:Ny-2] û[Ny-1:2*Ny-4] û[2*Ny-3:3*Ny-6]]
    u = zeros(Nz,Ny,Nx,3)
    for d=1:3, i=1:Nx, j=1:Ny-2, k=1:Nz
        u[k,j+1,i,d] = (û[j,d]*exp(im*(α*x[i]+β*z[k]))+conj(û[j,d])*exp(-im*(α*x[i]+β*z[k])))
    end
    return u
end

# Global Framework

using FFTW, SparseArrays

function uvec2ûvec_global(uvec,Nx,Ny,Nz)
    
    u = reshape(uvec,Nx,Ny,Nz,3)
    ûvec = permutedims(fft(u,(1,3))[:,2:Ny-1,:,:],(2,4,1,3))[:]
    
    return ûvec
end

function ûvec2uvec_global(ûvec,Nx,Ny,Nz)
    
    û = ûvec2û_global(ûvec,Nx,Ny,Nz)
    uvec = real.(ifft(û,(1,3)))[:]
    
    return uvec
end

function ûvec2û_global(ûvec,Nx,Ny,Nz)
    
    û = permutedims(reshape(ûvec,Ny-2,3,Nx,Nz),(3,1,4,2))
    wall = zeros(Nx,1,Nz,3)
    û = cat(wall,û,wall,dims=2)
    
    return û
end

function û2ûvec_global(û,Ny)
    
    ûvec = permutedims(û[:,2:Ny-1,:,:],(2,4,1,3))[:]
    
    return ûvec
end

function spectral_pad(û,Nx,Nz)

    nx, Ny, nz = size(û)[1:3]
    ûpad = zeros(Complex{Float64},Nx,Ny,Nz,3)

    ûpad[1:(nx+1)÷2, :, 1:(nz+1)÷2, :] = û[1:(nx+1)÷2, :, 1:(nz+1)÷2, :]
    ûpad[1:(nx+1)÷2, :, end-nz÷2+1:end, :] = û[1:(nx+1)÷2, :, (nz+1)÷2+1:end, :]
    ûpad[end-nx÷2+1:end, :, 1:(nz+1)÷2, :] = û[(nx+1)÷2+1:end, :, 1:(nz+1)÷2, :]
    ûpad[end-nx÷2+1:end, :, end-nz÷2+1:end, :] = û[(nx+1)÷2+1:end, :, (nz+1)÷2+1:end, :]

    return ûpad
end

function spectral_chop(û,nx,nz)

    Ny = size(û)[2]
    ûchop = zeros(Complex{Float64},nx,Ny,nz,3)

    ûchop[1:(nx+1)÷2, :, 1:(nz+1)÷2, :] = û[1:(nx+1)÷2, :, 1:(nz+1)÷2, :]
    ûchop[1:(nx+1)÷2, :, (nz+1)÷2+1:end, :] = û[1:(nx+1)÷2, :, end-nz÷2+1:end, :]
    ûchop[(nx+1)÷2+1:end, :, 1:(nz+1)÷2, :] = û[end-nx÷2+1:end, :, 1:(nz+1)÷2, :]
    ûchop[(nx+1)÷2+1:end, :, (nz+1)÷2+1:end, :] = û[end-nx÷2+1:end, :, end-nz÷2+1:end, :]

    return ûchop
end

# function pois_operator_global(Nx,Ny,Nz,R)

#     n = 2*Ny-4
#     A = spzeros(Complex{Float64},n*Nx*Nz,n*Nx*Nz)
#     B = spzeros(Complex{Float64},n*Nx*Nz,(3*Ny-6)*Nx*Nz)
#     C = spzeros(Complex{Float64},(3*Ny-6)*Nx*Nz,n*Nx*Nz)
#     Q = spzeros(Complex{Float64},n*Nx*Nz,n*Nx*Nz)


#     N = Ny
#     y = cheb_grid(N)
#     U = 1 .-y.^2 #base flow
#     D = cheb_diff(N)
#     D2 = D^2
#     S = Diagonal([0; 1 ./(1 .-y[2:N-1].^2); 0])
#     D4 = (Diagonal(1 .- y.^2)*D^4 - 8*Diagonal(y)*D^3 - 12*D^2)*S
#     W = Diagonal(cheb_weights(N))
#     Ũ = U[2:N-1]
#     D̃ = D[2:N-1,2:N-1]
#     D̃2 = D2[2:N-1,2:N-1]
#     D̃4 = D4[2:N-1,2:N-1]
#     W̃ = W[2:N-1,2:N-1]

#     α = ifftshift(-Nx÷2:(Nx+1)÷2-1)
#     β = ifftshift(-Nz÷2:(Nz+1)÷2-1)

#     for j=1:Nz, i=1:Nx
        
#         if α[i]==0 && β[j]==0
        
#             Amn = [D̃2/R zero(D̃);
#                 zero(D̃) D̃2/R]

#             Bmn = [I zero(D̃) zero(D̃);
#                 zero(D̃) zero(D̃) I]

#             Cmn = [I zero(D̃);
#                 zero(D̃) zero(D̃);
#                 zero(D̃) I]

#             Qmn = [W̃ zero(D̃);
#                 zero(D̃) W̃]/2
        
#         else
        
#             k = sqrt(α[i]^2+β[j]^2)

#             M = k^2*I-D̃2
#             M2 = k^4*I-2*k^2*D̃2+D̃4
#             invM = M\I
        
#             Aos = invM*(-im*α[i]*Diagonal(Ũ)*M-im*α[i]*Diagonal(D̃2*Ũ)-(1/R)*M2)
#             Asq = -im*α[i]*Diagonal(Ũ)-(1/R)*M
#             Ac = Array(-im*β[j]*Diagonal(D̃*Ũ))

#             Amn = [Aos zero(D̃);
#                 Ac Asq]
            
#             Bmn = [im*α[i]*invM*D̃ invM*k^2 im*β[j]*invM*D̃;
#                 im*β[j]*I zero(D̃) -im*α[i]*I]

#             Cmn = [im*α[i]/k^2*D̃ -im*β[j]/k^2*I;
#                  I zero(D̃);
#                  im*β[j]/k^2*D̃ im*α[i]/k^2*I]

#             Qmn = [(k^2*W̃ + D̃'*W̃*D̃) zero(D̃);
#                 zero(D̃) W̃]/(2*k^2)

#         end

#         m = (i-1)*n + (j-1)*n*Nx
#         m2 = (i-1)*(3*Ny-6) + (j-1)*(3*Ny-6)*Nx
#         A[m+1:m+n, m+1:m+n] = Amn
#         B[m+1:m+n, m2+1:m2+3*Ny-6] = Bmn
#         C[m2+1:m2+3*Ny-6, m+1:m+n] = Cmn
#         Q[m+1:m+n, m+1:m+n] = Qmn
    
#     end
        
#     return A,B,C,Q
    
# end

function pois_operator_global(Nx,Ny,Nz,R)

    n = 2*Ny-4
    A = spzeros(Complex{Float64},n*((Nx÷2)+1)*Nz,n*((Nx÷2)+1)*Nz)
    Aadj = spzeros(Complex{Float64},n*((Nx÷2)+1)*Nz,n*((Nx÷2)+1)*Nz)
    B = spzeros(Complex{Float64},n*((Nx÷2)+1)*Nz,(3*Ny-6)*((Nx÷2)+1)*Nz)
    C = spzeros(Complex{Float64},(3*Ny-6)*((Nx÷2)+1)*Nz,n*((Nx÷2)+1)*Nz)
    Q = spzeros(Complex{Float64},n*((Nx÷2)+1)*Nz,n*((Nx÷2)+1)*Nz)


    N = Ny
    y = cheb_grid(N)
    U = 1 .-y.^2 #base flow
    D = cheb_diff(N)
    D2 = D^2
    S = Diagonal([0; 1 ./(1 .-y[2:N-1].^2); 0])
    D4 = (Diagonal(1 .- y.^2)*D^4 - 8*Diagonal(y)*D^3 - 12*D^2)*S
    W = Diagonal(cheb_weights(N))
    Ũ = U[2:N-1]
    D̃ = D[2:N-1,2:N-1]
    D̃2 = D2[2:N-1,2:N-1]
    D̃4 = D4[2:N-1,2:N-1]
    W̃ = W[2:N-1,2:N-1]

    α = collect(-Nx÷2:(Nx+1)÷2-1)
    β = ifftshift(-Nz÷2:(Nz+1)÷2-1)

    for j=1:Nz, i=1:((Nx÷2)+1)
        
        if α[i]==0 && β[j]==0
        
            Amn = [D̃2/R zero(D̃);
                   zero(D̃) D̃2/R]

            Aadjmn = [D̃2'/R zero(D̃);
                      zero(D̃) D̃2'/R]

            Bmn = [I zero(D̃) zero(D̃);
                  zero(D̃) zero(D̃) I]

            Cmn = [I zero(D̃);
                zero(D̃) zero(D̃);
                zero(D̃) I]

            Qmn = [W̃ zero(D̃);
                   zero(D̃) W̃]/2
        
        else
        
            k = sqrt(α[i]^2+β[j]^2)

            M = k^2*I-D̃2
            M2 = k^4*I-2*k^2*D̃2+D̃4
            invM = M\I
        
            Aos = invM*(-im*α[i]*Diagonal(Ũ)*M-im*α[i]*Diagonal(D̃2*Ũ)-(1/R)*M2)
            Asq = -im*α[i]*Diagonal(Ũ)-(1/R)*M
            Ac = Array(-im*β[j]*Diagonal(D̃*Ũ))

            Amn = [Aos zero(D̃);
                   Ac Asq]

            Aosadj = invM*(im*α[i]*Diagonal(Ũ)*M-2*im*α[i]*Diagonal(D̃*Ũ)*D̃-(1/R)*M2)
            Asqadj = im*α[i]*Diagonal(Ũ)-(1/R)*M
            Acadj = invM*Array(im*β[j]*Diagonal(D̃*Ũ))

            Aadjmn = [Aosadj Acadj;
                      zero(D̃) Asqadj]
            
            Bmn = [im*α[i]*invM*D̃ invM*k^2 im*β[j]*invM*D̃;
                im*β[j]*I zero(D̃) -im*α[i]*I]

            Cmn = [im*α[i]/k^2*D̃ -im*β[j]/k^2*I;
                 I zero(D̃);
                 im*β[j]/k^2*D̃ im*α[i]/k^2*I]

            Qmn = 2*[(k^2*W̃ + D̃'*W̃*D̃) zero(D̃);
                zero(D̃) W̃]/(2*k^2)

        end

        m = (i-1)*n + (j-1)*n*((Nx÷2)+1)
        m2 = (i-1)*(3*Ny-6) + (j-1)*(3*Ny-6)*((Nx÷2)+1)
        A[m+1:m+n, m+1:m+n] = Amn
        Aadj[m+1:m+n, m+1:m+n] = Aadjmn
        B[m+1:m+n, m2+1:m2+3*Ny-6] = Bmn
        C[m2+1:m2+3*Ny-6, m+1:m+n] = Cmn
        Q[m+1:m+n, m+1:m+n] = Qmn
    
    end
        
    return A,Aadj,B,C,Q
    
end


function spectral_symmetry(ûvec,Nx,Ny,Nz)

    if length(ûvec) == (3*Ny-6)*((Nx÷2)+1)*Nz
    
    û = ûvec2û_global(ûvec,(Nx÷2)+1,Ny,Nz)
    û_sym = zeros(Complex{Float64},Nx,Ny,Nz,3)
    
    û_sym[1, :, :, :] = @view û[(Nx÷2)+1, :, :, :]
    û_sym[((Nx+1)÷2)+1:Nx, :, :, :] = @view û[1:(Nx÷2), :, :, :]
    û_sym[2:((Nx+1)÷2),:, :, :] = conj(@view û[(Nx÷2):-1:2-Nx%2, :, :, :])
        
    else
        
    û = ûvec2û_global(ûvec,Nx,Ny,Nz)
    û_sym = fftshift(û,1)[1:(Nx÷2)+1,:,:,:]
        
    end
       
    ûvec = û2ûvec_global(û_sym,Ny)
    return ûvec
end
