function [rho,u,p,E,grid] = make_grid

%%% Initialize memory %%%
%[DEFAULT] Setup Grid: (Boundary Grid):
grid.Nx = 100; % Only specified here
grid.xmin = 0;
grid.xmax = 1.0;
grid.time = 0;
grid.t_max = 0.2;
grid.Output_interval = 1000;
Nx = grid.Nx;

%[DEFAULT] Constants, updated in IC.m
grid.iter = 1;

%[DEFAULT] Grids, updated in IC.m
grid.dx = (grid.xmax - grid.xmin)/grid.Nx;
grid.time = 0;
grid.cfl = 0.98; %clf = udt/dx <= C_max
grid.dt = 0.98*grid.dx/100;
grid.NT = ceil(grid.t_max/grid.dt);

%Grid
grid.x = linspace(grid.xmin,grid.xmax,Nx);


%Setup
% JE2 (P1)
grid.density_left = 1.0; 
grid.velocity_left = 0.75;
grid.pressure_left = 1.0;
grid.density_right = 0.125;
grid.velocity_right = 0.0;
grid.pressure_right = 0.1;


%Quantities
rho = zeros(1,Nx);
p = zeros(1,Nx);
u = zeros(1,Nx);

%Right and left
grid.L = zeros(1,Nx);
grid.R = zeros(1,Nx);

%Setup IC:
i = 1;
while i <= Nx
    if i < Nx*0.33333333
        grid.L(i) = 1;
    else
        grid.R(i) = 1;
    end
    i = i + 1;
end
rho = rho + grid.L*grid.density_left;
rho = rho + grid.R*grid.density_right;
u = u + grid.L*grid.velocity_left;
u = u + grid.R*grid.velocity_right;
p = p + grid.L*grid.pressure_left;
p = p + grid.R*grid.pressure_right;


%Grid:
grid.gamma = 1.4;

%grid.specific_internal_energy (?)
E = p./( (grid.gamma -1).*rho) + 0.5*u.*u;

%Old option from wiki:
%grid.specific_internal_energy = mean(p./(rho.*(grid.gamma-1)));
%p = rho.*(grid.gamma-1)*grid.specific_internal_energy;
%E = p/(grid.gamma -1) + 0.5*rho.*u.*u;


%Lastly Print Stability / Stats
%fprintf("Stability Requires: (cdt/dx) we have C = %1.3f of max 1.0\n",grid.c*grid.dt/grid.dx);
fprintf("Grid: Nx: %d, NT: %d\n",grid.Nx,grid.NT);
fprintf("Grid-Spacing: dx: %g, dT: %g\n",grid.dx,grid.dt);


end