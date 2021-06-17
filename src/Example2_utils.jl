## --> Constants

R = 2000
E0 = 1e-5
Nx,Ny,Nz = [32,65,32]
x = fourier_grid(Nx)
y = cheb_grid(Ny)[Ny:-1:1]
z = fourier_grid(Nz)
dt = 0.5
uB = base_flow(Nx,Ny,Nz)
n = length(uB)
Q = energy_weight(reshape(uB,Nx,Ny,Nz,3))
F = factorize_weights(Q)
Fi = F\I

## --> Functions

function load_snapshots(dir,runs,t)
    m = length(t)
    nruns = length(runs)
    X = Array{Float64}(undef,n,(m-1)*nruns)
    Y = Array{Float64}(undef,n,(m-1)*nruns)
    for run in runs, j = 1:m
        u = (load_field(dir*"$(run)/u$(t[j])00.h5")[1] -uB)[:]/E0
        if j < m
            X[:,(run-1)*(m-1)+j] = F*u
        end
        if j > 1
            Y[:,(run-1)*(m-1)+j-1] = F*u
        end
    end
    return X,Y
end

function load_snapshots(dir,t)
    m = length(t)
    X = Array{Float64}(undef,n,(m-1))
    Y = Array{Float64}(undef,n,(m-1))
    for j = 1:m
        # u = (load_field(dir*"/u$(t[j])00.h5")[1] -uB)[:]/E0
        u = (load_field(dir*"/u$(t[j])00.h5")[1])[:]/E0
        if j < m
            X[:,j] = F*u
        end
        if j > 1
            Y[:,j-1] = F*u
        end
    end
    return X,Y
end

function save_mode(dir,phi,psi)
    dict = @strdict phi psi x y z Nx Ny Nz
    matwrite(dir,dict; compress=true)
end

function load_mode(dir)
    dict = matread(dir)
    ϕ = reshape(dict["phi"],Nx,Ny,Nz,3)
    ϕ = ϕ/energy(ϕ)
    ψ = reshape(dict["psi"],Nx,Ny,Nz,3)
    ψ = ψ/energy(ψ)
    return ϕ,ψ
end

function load_true_modes()
    dir = datadir("sims","channel_flow","resolvent_truth_w=0.0_mode1.mat")
    ϕtrue1,ψtrue1 = load_mode(dir)
    dir = datadir("sims","channel_flow","resolvent_truth_w=0.0_mode2.mat")
    ϕtrue2,ψtrue2 = load_mode(dir)
    dir = datadir("sims","channel_flow","resolvent_truth_w=0.0_mode3.mat")
    ϕtrue3,ψtrue3 = load_mode(dir)
    dir = datadir("sims","channel_flow","resolvent_truth_w=0.4_mode1.mat")
    ϕtrue4,ψtrue4 = load_mode(dir)
    ϕtrue = [ϕtrue1,ϕtrue2,ϕtrue3,ϕtrue4]
    ψtrue = [ψtrue1,ψtrue2,ψtrue3,ψtrue4]
    return ϕtrue,ψtrue
end

function align_mode(Φ,Ψ,ϕtrue,ψtrue,l,j)
    ϕdata = real.(Φ[:,l*(j-1)+1:l*j])
    ψdata = real.(Ψ[:,l*(j-1)+1:l*j])
    ϕdata = ϕdata/(ϕdata'*Q*ϕdata)*(ϕdata'*Q*(ϕtrue[:]))
    ψdata = ψdata/(ψdata'*Q*ψdata)*(ψdata'*Q*(ψtrue[:]))
    ϕdata = reshape(ϕdata,Nx,Ny,Nz,3)
    ψdata = reshape(ψdata,Nx,Ny,Nz,3)
    return ϕdata,ψdata
end

function get_modes_aligned(λ̃,Ṽ,ϕtrue,ψtrue)
    ϕtrue1,ϕtrue2,ϕtrue3,ϕtrue4 = ϕtrue
    ψtrue1,ψtrue2,ψtrue3,ψtrue4 = ψtrue
    ω = 0.0
    Ψ,σ,Φ = opt_forcing(λ̃,Ṽ,Q,ω)
    ϕdata1,ψdata1 = align_mode(Φ,Ψ,ϕtrue1,ψtrue1,2,1)
    ϕdata2,ψdata2 = align_mode(Φ,Ψ,ϕtrue2,ψtrue2,2,2)
    ϕdata3,ψdata3 = align_mode(Φ,Ψ,ϕtrue3,ψtrue3,2,3)
    σ1 = σ[1]
    σ2 = σ[3]
    σ3 = σ[5]
    ω = 0.4
    Ψ,σ,Φ = opt_forcing(λ̃,Ṽ,Q,ω)
    ϕdata4,ψdata4 = align_mode(Φ,Ψ,ϕtrue4,ψtrue4,4,1)
    σ4 = σ[1]
    ϕdata = [ϕdata1,ϕdata2,ϕdata3,ϕdata4]
    ψdata = [ψdata1,ψdata2,ψdata3,ψdata4]
    σdata = [σ1,σ2,σ3,σ4]
    return ϕdata,ψdata,σdata
end

function mode_error(ϕdata,ψdata,ϕtrue,ψtrue)
    ϕtrue1,ϕtrue2,ϕtrue3,ϕtrue4 = ϕtrue
    ψtrue1,ψtrue2,ψtrue3,ψtrue4 = ψtrue
    ϕdata1,ϕdata2,ϕdata3,ϕdata4 = ϕdata
    ψdata1,ψdata2,ψdata3,ψdata4 = ψdata
    eψ1 = energy(ψtrue1-ψdata1)
    eϕ1 = energy(ϕtrue1-ϕdata1)
    eψ2 = energy(ψtrue2-ψdata2)
    eϕ2 = energy(ϕtrue2-ϕdata2)
    eψ3 = energy(ψtrue3-ψdata3)
    eϕ3 = energy(ϕtrue3-ϕdata3)
    eψ4 = energy(ψtrue4-ψdata4)
    eϕ4 = energy(ϕtrue4-ϕdata4)
    eϕ = [eϕ1,eϕ2,eϕ3,eϕ4]
    eψ = [eψ1,eψ2,eψ3,eψ4]
    return eϕ,eψ
end

using ApproxFun

function u2û(u)
    nβ = Nz-2
    ~,~,~,~,~,u_yz = u_slices(u,1)
    ~,~,~,~,~,v_yz = u_slices(u,2)
    ~,~,~,~,~,w_yz = u_slices(u,3)
    S = Laurent()
    û = zeros(Complex{Float64},3*Ny-6,nβ)
    for j=1:Ny-2
        û[j,:] = coefficients(Fun(S,ApproxFun.transform(S,u_yz[:,j+1])))[2:nβ+1]
        û[Ny-2+j,:] = coefficients(Fun(S,ApproxFun.transform(S,v_yz[:,j+1])))[2:nβ+1]
        û[2*Ny-4+j,:] = coefficients(Fun(S,ApproxFun.transform(S,w_yz[:,j+1])))[2:nβ+1]
    end
    return û
end

function q̂_response(f̂,ω,p)
    N = div(length(f̂),2)+2
    A,Q = pois_operator(N,p)
    q̂ = inv(-im*ω*I-A)*f̂
    return q̂
end

function u_response(f)
    f̂u = u2û(f)
    nβ = Nz-2
    u = zero(uB)
    β = zeros(nβ)
    β[1:2:end] = -collect(1:nβ/2)
    β[2:2:end] = collect(1:nβ/2)
    f̂q = zeros(Complex{Float64},2*Ny-4,nβ)
    q̂ = zeros(Complex{Float64},2*Ny-4,nβ)
    for k=1:nβ
        f̂q[:,k] = u2q(f̂u[:,k],0.0,β[k])
        q̂[:,k] = q̂_response(f̂q[:,k],0,[0,β[k],R])
        u += q̂2u(q̂[:,k],0,β[k],Nx,Nz) / 2
    end
    return u
end
;
