using DrWatson
@quickactivate "Data-driven resolvent"
using LinearAlgebra, ApproxFun, DifferentialEquations
include(srcdir("NonModalStabilityTools.jl"))
include(srcdir("SpectralDiff.jl"))

function gl_operator(N::Integer, p)
    μ₀, μ₂, ν, γ = p
    cᵤ = imag(ν)/2
    χ = (-μ₂/(2γ))^0.25
    b = real.(χ)
    x = her_grid(N,b)
    μ = μ₀-cᵤ^2 .+μ₂*x.^2/2
    ∂ₓ = her_diff(N,b,1)
    ∂ₓₓ = her_diff(N,b,2)
    A = -ν*∂ₓ+γ*∂ₓₓ+Diagonal(μ)
    Q = Diagonal(her_weights(N,b))
    return A,Q
end

function gl_subcritical_setup()
    N = 220
    μ₀ = 0.23
    μ₂ = -0.01
    ν = 2.0+0.4*im
    γ = 1.0 - 1.0*im
    p = [μ₀, μ₂, ν, γ]
    χ = (-μ₂/(2γ))^0.25
    b = real.(χ)
    x = her_grid(N,b)
    A,Q = gl_operator(N,p)
    return A,Q,x,p
end

function gl_ic_basis(howmany)
    A,Q,x,p = gl_subcritical_setup()
    N = length(x)
    μ₀, μ₂, ν, γ = p
    χ = (-μ₂/(2γ))^0.25
    b = real.(χ)
    S = GaussWeight(Hermite(b^2),b^2/2)
    q0 = hcat([(Fun(S,I(howmany)[i,:]).(x) .+0*im) for i=1:howmany]...)
    q0 = normalize_basis(q0,Q)
    return q0
end

function gl_run_transient(howmany;verbose=true)
    A,Q,x,p = gl_subcritical_setup()
    dt = 0.5
    T = 50.0
    t = collect(0.0:dt:T)
    m = length(t)
    tol=1e-16
    gl_ode(dq,q,p,t) = mul!(dq,A,q)
    q = []
    if howmany==0
        verbose ? println("solving for optimal disturbance") : nothing
        q0 = opt_disturbance(A,Q,8.25)[3][:,1]
        prob = ODEProblem(gl_ode,q0,(0,T))
        push!(q,Array(solve(prob,saveat=dt, reltol=tol, abstol=tol)))
    else
        q0 = gl_ic_basis(howmany)
        for i=1:howmany
            verbose ? println("running transient $i")  : nothing
            prob = ODEProblem(gl_ode,q0[:,i],(0,T))
            push!(q,Array(solve(prob,saveat=dt, reltol=tol, abstol=tol)))
        end
    end
    return t,x,q
end
