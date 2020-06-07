function dy = dyGrowth(t,y,T,D,M,L,ka,rho,gamma,Ntot,dp0)
% Differential equation system to be solved by ode15s solver. 
%
% INPUT
% t - time (s)
% y - y(1) particle diameter (m) and y(2) partial pressure of water (Pa)
% T - Temperature (K)
% D - Diffusion coefficient of water, (m^2/s)
% M - Molar mass of water (kg/mol)
% L - Heat of evaporation of water (J/kg)
% ka - Thermal conductivity of air (W/(m*K))
% rho - Density of water (kg/m^3)
% gamma - Surface tension of water (N/m)
% Ntot - Number density of particles (#/m^3)
% dp0 - particle diameter in the beginning (m)
%
% OUTPUT
% dy - rates of changes. dy(1) particle diameter (m) and 
% dy(2) partial pressure of water (Pa) 

dy = zeros(2,1);
Td = droplet_temp(T, D, M, L, ka, y(2),y(1),rho,gamma);
dy(1) = dp_growth_rate( y(1), T, rho, gamma, M, D, y(2), Td);
dy(2) = p_growth_rate( T, rho, Ntot, y(1), dy(1), M );

% If particle would like to evaporate, but it is at the original size,
% the diameter is not changed smaller. We assume non-evaporating nucleus.

if dy(1) < 0 && y(1) <= dp0
    dy(1) = 0;
    dy(2) = 0;
end

end

