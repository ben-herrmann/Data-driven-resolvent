using DrWatson
@quickactivate "Data-driven resolvent"
using LinearAlgebra, ApproxFun

# fourier

function fourier_grid(N::Integer)
    x = collect(0:N-1)*2π/N
    return x
end

function fourier_diff(N::Integer,degree::Integer=1)
    Δx=2π/N
    if degree%2≠0
        D1 = zeros(N,N)
        for i=1:N, j=1:N
            if i≠j
                D1[i,j] = 0.5*(-1)^(i+j)*cot(0.5*(i-j)*Δx)
            end
        end
    elseif degree>=2
        D2 = fill(-π^2/(3*Δx^2)-1/6,(N,N))
        for i=1:N, j=1:N
            if i≠j
                D2[i,j] = -0.5*(-1)^(i+j)*csc(0.5*(i-j)*Δx)^2
            end
        end
    end
    if degree==1
        D=D1
    elseif degree%2==0
        D=D2^(degree/2)
    else
        D = D2^div(degree,2)*D1^(degree%2)
    end
    return D
end

function fourier_weights(N::Integer)
    w = 2π/N*ones(N)
    return w
end

function fourier_interp(v::Array,xq)
    coeffs = ApproxFun.transform(Fourier(),v)
    vq = Fun(Fourier(),coeffs).(xq)
    return vq
end

# chebyshev

function cheb_grid(N::Integer)
    x = cos.(π.*(0:N-1)./(N-1))
    return x
end

function cheb_diff(N::Integer)
    x = cheb_grid(N)
    c = [2; ones(N-2); 2].*(-1).^(0:N-1)
    X = repeat(x,1,N)
    dX = X-X'
    D = (c*(1 ./c)')./(dX+I)
    D = D-Diagonal(vec(sum(D,dims=2)))
    return D
end

function cheb_weights(N::Integer)
    θ = π*(0:N-1)/(N-1)
    x = cos.(θ)
    w = zeros(N)
    v = ones(N-2)
    if N%2 ≠ 0
        w[1] = 1/((N-1)^2-1)
        w[N] = w[1]
        for k=1:(N-1)/2-1
            v = v.-2*cos.(2*k*θ[2:N-1])/(4*k^2-1)
        end
        v = v.-cos.((N-1)*θ[2:N-1])/((N-1)^2-1)
    else
        w[1] = 1/(N-1)^2
        w[N] = w[1]
        for k=1:(N-2)/2
            v = v.-2*cos.(2*k*θ[2:N-1])/(4*k^2-1)
        end
    end
    w[2:N-1] = 2*v/(N-1)
    return w
end

function cheb_interp(v::Array,xq,m=50)
    S = Chebyshev()
    n = length(v)
    x = cheb_grid(n)
    basis = Array{Float64}(undef,n,m)
    for k = 1:m
        basis[:,k] = Fun(S,[zeros(k-1);1]).(x)
    end
    vq = Fun(S,basis\v).(xq)
    return vq
end

# gauss-weghted hermite

function her_grid(N::Integer,b::Real=1.0)
    #  J.A.C. Weideman, S.C. Reddy 1998.
    J = SymTridiagonal(zeros(N),sqrt.(1:N-1))
    r = eigvals(J)/sqrt(2)
    x = r/b
    return x
end

function her_diff(N::Integer,b::Real=1.0,degree::Integer=1.0)
    #  J.A.C. Weideman, S.C. Reddy 1998.
    x = her_grid(N)
    α = exp.(-x.^2/2)
    β = Array{typeof(x[1])}(undef,degree+1,N)
    β[1,:] = ones(length(x))
    β[2,:] = -x
    for l=3:degree+1
        β[l,:] = -x.*β[l-1,:].-(l-2)*β[l-2,:]
    end
    β = β[2:end,:]
    XX = repeat(x,1,N)
    DX = XX-XX'
    DX[I(N)].=1.0
    c = α.*prod(DX,dims=2)
    C = repeat(c,1,N)
    C = C./C'
    Z = 1 ./DX
    Z[I(N)].=0.0
    X = Z'
    X = reshape(X[.~I(N)],N-1,N)
    Y = ones(N-1,N)
    D = I(N)
    for d=1:degree
        Y = cumsum([β[d,:]'; d*Y[1:N-1,:].*X],dims=1)
        D = d*Z.*(C.*repeat(diag(D),1,N)-D)
        D[I(N)]=Y[N,:]
    end
    D = b^degree*D
    return D
end

function her_weights(N::Integer,b::Real=1.0)
    x = her_grid(N,b)
    w = ([diff(x);0]+[0;diff(x)])/2
    return w
end

# 2D-domains

## grids

function meshgrid(x, y)
    Nx = length(x)
    Ny = length(y)
    X = [x[i] for i=1:Nx,j=1:Ny]
    Y = [y[j] for i=1:Nx,j=1:Ny]
    return X, Y
end

function fourier_grid2D(Nx,Ny)
    x = fourier_grid(Nx)
    y = fourier_grid(Ny)
    X,Y = meshgrid(x,y)
    return X,Y
end

function fourierXcheb_grid(Nx,Ny)
    x = fourier_grid(Nx)
    y = cheb_grid(Ny)
    X,Y = meshgrid(x,y)
    return X,Y
end

function cheb_grid2D(Nx,Ny)
    x = cheb_grid(Nx)
    y = cheb_grid(Ny)
    X,Y = meshgrid(x,y)
    return X,Y
end

## interpolation

function fourier_interp2D(v::Array,x::Array,y::Array)
    S = Fourier()*Fourier()
    coeffs = ApproxFun.transform(S,v[:])
    vq = Fun(S,coeffs).(x,y)
    return vq
end

function fourier_interp2D(v::Array,Nx::Int,Ny::Int)
    x,y = fourier_grid2D(Nx,Ny)
    vq = fourier_interp2D(v,x,y)
    return vq
end

function fourierXcheb_interp(v::Array,xq::Array,yq::Array,m=50)
    S = Fourier()*Chebyshev()
    x, y = fourierXcheb_grid(size(v)...)
    n = length(v)
    basis = Array{Float64}(undef,n,m)
    for k = 1:m
        basis[:,k] = (Fun(S,[zeros(k-1);1]).(x,y))[:]
    end
    vq = Fun(S,basis\(v[:])).(xq,yq)
    return vq
end

function fourierXcheb_interp(v::Array,Nxq::Int,Nyq::Int,m=50)
    xq, yq = fourierXcheb_grid(Nxq,Nyq)
    vq = fourierXcheb_interp(v,xq,yq,m)
end

function cheb_interp2D(v::Array,xq::Array,yq::Array,m=50)
    S = Chebyshev()*Chebyshev()
    x, y = cheb_grid2D(size(v)...)
    n = length(v)
    basis = Array{Float64}(undef,n,m)
    for k = 1:m
        basis[:,k] = (Fun(S,[zeros(k-1);1]).(x,y))[:]
    end
    vq = Fun(S,basis\(v[:])).(xq,yq)
    return vq
end

function cheb_interp2D(v::Array,Nxq::Int,Nyq::Int,m=50)
    xq, yq = cheb_grid2D(Nxq,Nyq)
    vq = cheb_interp2D(v,xq,yq,m)
end