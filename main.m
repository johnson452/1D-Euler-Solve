%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: Grant Johnson
%Date: 5/30/2022
%Yee algorithm *1D*
%Non-relativistic Euler, verify MUSCL scheme

%Notes:
%-1D
% Fluid Scheme: https://ammar-hakim.org/sj/hancock-muscl.html
% Good example: https://en.wikipedia.org/wiki/MUSCL_scheme
% https://www.cambridge.org/core/services/aop-cambridge-core/content/view/8F5CD408E7073099BFDE1E409C1E79AB/S0022112077001463a.pdf/a-numerical-study-of-a-converging-cylindrical-shock.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%Setup Initial shock-tube problem 
[rho,u,p,E,grid] = make_grid();

%Make the diagnostic Figure
figure('units','normalized','outerposition',[0 0 0.5 0.5])

%%% Time loop %%%
while(grid.time < grid.t_max)
    
    %Call i/o and diagnostics
    diagnostics(rho,u,p,E,grid);
    
    %Update the gridtime
    grid.time = grid.time + grid.dt;
    
    %Update the iterator
    grid.iter = grid.iter + 1;
 
    %Updater - updates all quantities simultaneosly
    % n -> n + 1 all quantities
    [rho,u,p,E,grid] = push_all(rho,u,p,E,grid);

     %BC - All outflow (copy)
    [rho,u,p,E,grid] = BC(rho,u,p,E,grid);

end
%%% End Time Loop %%%
%%% End main %%%