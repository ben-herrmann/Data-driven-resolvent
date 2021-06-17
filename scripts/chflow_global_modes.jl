using DrWatson
@quickactivate "Data-driven resolvent"

include(srcdir("PoiseuilleFlow.jl"))

Nx = 128
Nz = 128
data = wload(datadir("sims","channel_flow","spectrum.bson"))
@unpack λ,V̂,Ŵ,α,β = data
r = 10
V = hcat([q̂2u(V̂[:,i],α[i],β[i],Nx,Nz)[:] for i=1:r]...)
W = hcat([q̂2u(Ŵ[:,i],α[i],β[i],Nx,Nz)[:] for i=1:r]...)
data = @dict(λ,V,W,α,β)
wsave(datadir("sims","channel_flow","global_modes.bson"),data)