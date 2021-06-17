using DrWatson
@quickactivate "Data-driven resolvent"
using LinearAlgebra

function adj(A,Q)
    Aadj = Q\A'*Q
    return Aadj
end

function normalize_basis(V,Q)
    for i=1:size(V,2)
        V[:,i] = V[:,i]/sqrt(V[:,i]'*Q*V[:,i])
    end
    return V
end

function eigen_dual(A,Q,log_sort::Bool=false)
    Aadj = adj(A,Q)
    if log_sort
        λ, V = eigen(A, sortby=x->imag(log(x)))
        λ̄, W = eigen(Aadj, sortby=x->-imag(log(x)))
        p = sortperm(-real(log.(λ)))
        p̄ = sortperm(-real(log.(λ̄)))
    else
        λ, V = eigen(A, sortby=x->imag(x))
        λ̄, W = eigen(Aadj, sortby=x->-imag(x))
        p = sortperm(-real(λ))
        p̄ = sortperm(-real(λ̄))
    end
    V = V[:,p]
    λ = λ[p]
    W = W[:,p̄]
    λ̄ = λ̄[p̄]
    V = normalize_basis(V,Q)
    W = normalize_basis(W,Q)
    for i=1:size(V,2)
        V[:,i] = V[:,i]/(W[:,i]'*Q*V[:,i])
    end
    return λ, V, W
end

function factorize_weights(Q)
    if typeof(Q)<:Diagonal
        F = sqrt.(Q)
    else
        F = cholesky(Hermitian(Q)).U
    end
    return F
end

function opt_disturbance(A,Q,t)
    F = factorize_weights(Q)
    Fi = F\I
    Ψ,s,Φ = svd(F*exp(A*t)*Fi,full=false)
    Ψ = Fi*Ψ
    Φ = Fi*Φ
    return Ψ,s,Φ
end

function opt_disturbance(λ,V,Q,t)
    Ã = Diagonal(λ)
    Q̃ = Array(V'*Q*V)
    Ψ,σ,Φ = opt_disturbance(Ã,Q̃,t)
    Ψ = V*Ψ
    Φ = V*Φ
    return Ψ,σ,Φ
end

function opt_growth(A,Q,tspan,m)
    F = factorize_weights(Q)
    Fi = F\I
    G = []
    for t in tspan
        push!(G,(svdvals(F*exp(A*t)*Fi)[1:m]).^2)
    end
    G = hcat(G...)
    return G
end

function opt_growth(λ,V,Q,ωspan,m)
    Ã = Diagonal(λ)
    Q̃ = Array(V'*Q*V)
    G = opt_growth(Ã,Q̃,ωspan,m)
    return G
end

function opt_forcing(A,Q,ω)
    F = factorize_weights(Q)
    Fi = F\I
    Ψ,σ,Φ = svd(F*inv(-im*ω*I-A)*Fi,full=false)
    Ψ = Fi*Ψ
    Φ = Fi*Φ
    return Ψ,σ,Φ
end

function opt_forcing(λ,V,Q,ω)
    Ã = Diagonal(λ)
    Q̃ = Array(V'*Q*V)
    Ψ,σ,Φ = opt_forcing(Ã,Q̃,ω)
    Ψ = V*Ψ
    Φ = V*Φ
    return Ψ,σ,Φ
end

function opt_forcing(C,A,B,Qy,Qu,ω)
    Fy = factorize_weights(Qy)
    Fyi = Fy\I
    Fu = factorize_weights(Qu)
    Fui = Fu\I
    Ψ,σ,Φ = svd(Fy*C*inv(-im*ω*I-A)*B*Fui,full=false)
    Ψ = Fyi*Ψ
    Φ = Fui*Φ
    return Ψ,σ,Φ
end

function opt_forcing(C,λ,V,W,B,Qy,Q,Qu,ω)
    C̃ = C*V
    Ã = Diagonal(λ)
    B̃ = W'*Q*B
    Ψ,σ,Φ = opt_forcing(C̃,Ã,B̃,Qy,Qu,ω)
    return Ψ,σ,Φ
end

function opt_gain(A,Q,ωspan,m)
    F = factorize_weights(Q)
    Fi = F\I
    R = []
    for ω in ωspan
        push!(R,(svdvals(F*inv(-im*ω*I-A)*Fi)[1:m]).^2)
    end
    R = hcat(R...)
    return R
end

function opt_gain(λ,V,Q,ωspan,m)
    Ã = Diagonal(λ)
    Q̃ = Array(V'*Q*V)
    R = opt_gain(Ã,Q̃,ωspan,m)
    return R
end

function opt_gain(C,A,B,Qy,Qu,ωspan,m)
    Fy = factorize_weights(Qy)
    Fu = factorize_weights(Qu)
    Fui = Fu\I
    R = []
    for ω in ωspan
        push!(R,(svdvals(Fy*C*inv(-im*ω*I-A)*B*Fui)[1:m]).^2)
    end
    R = hcat(R...)
    return R
end

function opt_gain(C,λ,V,W,B,Qy,Q,Qu,ωspan,m)
    C̃ = C*V
    Ã = Diagonal(λ)
    B̃ = W'*Q*B
    R = opt_gain(C̃,Ã,B̃,Qy,Qu,ωspan,m)
    return R
end

function DMD(X,Y,dt,r,tol=1e-6)
    nt = size(X,2)
    U, s, V = svd(X, full=false)
    r = minimum((r, size(U,2)))
    Uᵣ = @view U[:,1:r]
    Sᵣ = @view Diagonal(s)[1:r,1:r]
    Vᵣ = @view V[:,1:r]
    Ã = (Uᵣ'*Y)*(Vᵣ/Sᵣ)
    ρ, W, Wadj = eigen_dual(Ã, Array(I(r)), true)
    Ψ = Y*(Vᵣ/Sᵣ*W)
    Φ = Uᵣ*Wadj
    @inbounds for i=1:r
        Ψ[:,i] = Ψ[:,i]/sqrt(Ψ[:,i]'*Ψ[:,i])
        Φ[:,i] = Φ[:,i]/sqrt(Φ[:,i]'*Φ[:,i])
        Ψ[:,i] = Ψ[:,i]/(Φ[:,i]'*Ψ[:,i])
    end
    b = Ψ\(X[:,1])
    large = abs.(b).>tol*maximum(abs.(b))
    Ψ = Ψ[:,large]
    Φ = Φ[:,large]
    ρ = ρ[large]
    λ = log.(ρ)/dt
    b = b[large]
    return λ, Ψ, Φ, b
end

function DMD_from_POD(X,Y,dt,U,s,V,r,tol=1e-16)
    nt = size(X,2)
    Uᵣ = @view U[:,1:r]
    Sᵣ = @view Diagonal(s)[1:r,1:r]
    Vᵣ = @view V[:,1:r]
    Ã = (Uᵣ'*Y)*(Vᵣ/Sᵣ)
    ρ, W, Wadj = eigen_dual(Ã, Array(I(r)), true)
    Ψ = Y*(Vᵣ/Sᵣ*W)
    Φ = Uᵣ*Wadj
    @inbounds for i=1:r
        Ψ[:,i] = Ψ[:,i]/sqrt(Ψ[:,i]'*Ψ[:,i])
        Φ[:,i] = Φ[:,i]/sqrt(Φ[:,i]'*Φ[:,i])
        Ψ[:,i] = Ψ[:,i]/(Φ[:,i]'*Ψ[:,i])
    end
    b = Ψ\(X[:,1])
    large = abs.(b).>tol*maximum(abs.(b))
    Ψ = Ψ[:,large]
    Φ = Φ[:,large]
    ρ = ρ[large]
    λ = log.(ρ)/dt
    b = b[large]
    return λ, Ψ, Φ, b
end

## For sparse matrices

using KrylovKit

function opt_gain_sparse(λ,V,Q,ωspan)
    Ã = Diagonal(λ)
    Q̃ = Array(V'*Q*V)
    F̃ = factorize_weights(Q̃)
    F̃i = sparse(F̃\I)
    F̃ = sparse(F̃)
    R = []
    for ω in ωspan
        push!(R,(svdsolve(F̃*inv(-im*ω*I-Ã)*F̃i,1,:LR)[1][1]).^2)
    end
    R = hcat(R...)
    return R
end

function opt_gain_sparse(C,λ,V,W,B,Qy,Q,Qu,ωspan)
    C̃ = C*V
    Ã = Diagonal(λ)
    B̃ = W'*Q*B
    Fy = factorize_weights(Array(Qy))
    Fyi = sparse(Fy\I)
    Fy = sparse(Fy)
    Fu = factorize_weights(Array(Qu))
    Fui = sparse(Fu\I)
    Fu = sparse(Fu)
    R = []
    for ω in ωspan
        push!(R,(svdsolve(Fy*C̃*inv(-im*ω*I-Ã)*B̃*Fui,1,:LR)[1][1]).^2)
    end
    R = hcat(R...)
    return R
end

function opt_forcing_sparse(λ,V,Q,ω)
    Ã = Diagonal(λ)
    Q̃ = Array(V'*Q*V)
    F̃ = factorize_weights(Q̃)
    F̃i = sparse(F̃\I)
    F̃ = sparse(F̃)
    σ, Ψ, Φ, ~ = svdsolve(F̃*inv(-im*ω*I-Ã)*F̃i)
    Φ = hcat(Φ...)
    Ψ = hcat(Ψ...)
    Ψ = V*F̃i*Ψ
    Φ = V*F̃i*Φ
    return Ψ,σ,Φ
end

function opt_forcing_sparse(C,λ,V,W,B,Qy,Q,Qu,ω)
    C̃ = C*V
    Ã = Diagonal(λ)
    B̃ = W'*Q*B
    Fy = factorize_weights(Array(Qy))
    Fyi = sparse(Fy\I)
    Fy = sparse(Fy)
    Fu = factorize_weights(Array(Qu))
    Fui = sparse(Fu\I)
    Fu = sparse(Fu)
    σ, Ψ, Φ, ~ = svdsolve(Fy*C̃*inv(-im*ω*I-Ã)*B̃*Fui)
    Φ = hcat(Φ...)
    Ψ = hcat(Ψ...)
    Ψ = Fyi*Ψ
    Φ = Fui*Φ
    return Ψ,σ,Φ
end
