using DrWatson
@quickactivate "Data-driven resolvent"

include(srcdir("SpectralDiff.jl"))
using HDF5

function base_flow(Nx,Ny,Nz)
    Nx = Int(Nx)
    Ny = Int(Ny)
    Nz = Int(Nz)
    uB = zeros(Nx,Ny,Nz,3)
    y = cheb_grid(Ny)
    for ix=1:Nx, iz=1:Nz
        uB[ix,:,iz,1] = 1 .-y.^2
    end
    return uB
end

function grid_coords(u)
    # Lx, Ly, Lz = [2π, 2, 2π]
    Nx, Ny, Nz = size(u)[1:3]
    x = fourier_grid(Nx)
    y = cheb_grid(Ny)[Ny:-1:1]
    z = fourier_grid(Nz)
    return x,y,z
end

function save_field(dir,u,attr)
    x,y,z = grid_coords(u)
    file = h5open(dir,"w")
    h5write(dir,"data/u",u)
    h5write(dir,"geom/x",x)
    h5write(dir,"geom/y",y)
    h5write(dir,"geom/z",z)
    h5writeattr(dir,"/",attr)
    close(file)
    return nothing
end

function load_field(dir)
    file = h5open(dir,"r")
    u = h5read(dir,"data/u")
    x = h5read(dir,"geom/x")
    y = h5read(dir,"geom/y")
    z = h5read(dir,"geom/z")
    attr = h5readattr(dir,"/")
    close(file)
    return u,x,y,z,attr
end

function energy_weight(u)
    Nx, Ny, Nz, Nd = size(u)
    Qx = Diagonal(fourier_weights(Nx))
    Qy = Diagonal(cheb_weights(Ny))
    Qz = Diagonal(fourier_weights(Nz))
    Q = kron(I(Nd),Qz,Qy,Qx)/(8π^2)
    return Q
end

function energy(u)
    Q = energy_weight(u)
    e = sqrt((u[:])'*Q*(u[:]))
    return e
end

function case_dir(case)
    E0,Nx,Ny,Nz = case
    Nx = Int(Nx)
    Ny = Int(Ny)
    dir = "E0=$(E0)_Nx=$(Nx)_Ny=$(Ny)"
    return dir
 end

function save_initial_condition(dir,case)
    E0,Nx,Ny,Nz = case
    attr = Dict("nu" => 0.0,"b" => 1.0,"a" => -1.0,
        "Lx" => 6.283185307179586,"Lz" => 6.283185307179586,
        "Nz" => Nz,"Nx" => Nx,"Ny" => Ny,"Nd" => 3,
        "Nypad" => Ny,"Nxpad" => Nx,"Nzpad" => Nz)
    up = localized_impulse(E0,Nx,Ny,Nz)
    uB = base_flow(Nx,Ny,Nz)
    subdir = case_dir(case)
    save_field(dir*subdir*"/uB.h5",uB,attr)
    save_field(dir*subdir*"/u0.h5",uB+up,attr)
    return nothing
end

function u_slices(u::Array,d::Int,Nxq::Int,Nyq::Int,Nzq::Int)
    u = u[:,:,:,d]
    Nz,Ny,Nx = size(u)
    u_zx = fourier_interp2D(Array(u[:,div(Ny+1,2),:]'),Nxq,Nzq)
    u_yx = Array(fourierXcheb_interp(Array(u[div(Nz,2),:,:]'),Nxq,Nyq)')
    u_zy = fourierXcheb_interp(u[:,:,div(Nx,2)],Nzq,Nyq)
    xq = fourier_grid(Nxq)
    yq = cheb_grid(Nyq)
    zq = fourier_grid(Nzq)
    return xq,yq,zq,u_zx,u_yx,u_zy
end

function u_slices(u::Array,d::Int)
    u = u[:,:,:,d]
    Nz,Ny,Nx = size(u)
    u_zx = Array(u[:,div(Ny+1,2),:]')
    u_xy = u[div(Nz,2),:,:]
    u_zy = Array(u[:,:,div(Nx,2)]')
    x = fourier_grid(Nx)
    y = cheb_grid(Ny)
    z = fourier_grid(Nz)
    return x,y,z,u_zx,u_xy,u_zy
end

function u_slices(u::Array,d::Int)
    u = u[:,:,:,d]
    Nz,Ny,Nx = size(u)
    u_zx = Array(u[:,div(Ny+1,2),:]')
    u_yx = u[div(Nz,2),:,:]
    u_zy = u[:,:,div(Nx,2)]
    x = fourier_grid(Nx)
    y = cheb_grid(Ny)
    z = fourier_grid(Nz)
    return x,y,z,u_zx,u_yx,u_zy
end
