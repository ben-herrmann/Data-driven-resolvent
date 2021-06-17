function U = get_component(u,l,Nx,Ny,Nz)
u = reshape(u,Nz,Ny,Nx,3);
U = zeros(Nx,Nz,Ny);
for i=1:Nx
    for j=1:Ny
        for k=1:Nz
            U(i,k,j) = u(i,j,k,l);
        end
    end
end
U = U/max(abs(u(:)));
end