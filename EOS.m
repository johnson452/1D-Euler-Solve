function [p] = EOS(rho,u,E,grid,type)

%Depending on the prescription return the EOS

if type == "calorically_ideal"
    p = (grid.gamma - 1)*rho.*(E - (1/2).*u.*u);
else
    p = rho.*(grid.gamma-1)*grid.specific_internal_energy;
end


end