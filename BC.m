function [rho,u,p,E,grid] = BC(rho,u,p,E,grid)

%Apply outflow BC on all quantities
Nx = grid.Nx;
rho(1) = rho(2);
u(1) = u(2);
p(1) = p(2);
E(1) = E(2);
rho(Nx) = rho(Nx-1);
u(Nx) = u(Nx-1);
p(Nx) = p(Nx-1);
E(Nx) = E(Nx-1);

end