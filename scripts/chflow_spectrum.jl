using DrWatson
@quickactivate "Data-driven resolvent"

include(srcdir("PoiseuilleFlow.jl"))

Ny = 101
αtmp = collect(-7:7)
βtmp = collect(-7:7)
R = 2000
allparams = @dict αtmp βtmp R
dicts = dict_list(allparams)
λ = []
V̂ = []
Ŵ = []
α = []
β = []
for dict in dicts
    @unpack αtmp,βtmp,R = dict
    if αtmp≠0 || βtmp≠0
        Atmp,Qtmp = pois_operator(Ny,[αtmp,βtmp,R])
        λtmp,Vtmp,Wtmp = eigen_dual(Atmp,Qtmp)
        push!(λ,λtmp)
        push!(V̂,Vtmp)
        push!(Ŵ,Wtmp)
        push!(α,fill(αtmp,length(λtmp)))
        push!(β,fill(βtmp,length(λtmp)))
    end
end
λ = vcat(λ...)
V̂ = hcat(V̂...)
Ŵ = hcat(Ŵ...)
α = vcat(α...).+0.0
β = vcat(β...).+0.0
perm = sortperm(-real(λ))
λ = λ[perm]
V̂ = V̂[:,perm]
Ŵ = Ŵ[:,perm]
α = α[perm]
β = β[perm]
data = @dict(λ,V̂,Ŵ,α,β)
wsave(datadir("sims","channel_flow","spectrum.bson"),data)