function [ Td ] = droplet_temp(T, D, M, L, ka, p,dp,rho,gamma  )
% Calculates surface temperature of particle during the condensation growth
%
% INPUT
% T - Temperature (K)
% D - Diffusion coefficient of water, (m^2/s)
% M - Molar mass of water (kg/mol)
% L - Heat of evaporation of water (J/kg)
% ka - Thermal conductivity of air (W/(m*K))
% p - Surrounding partial pressure of water (Pa)
% dp - particle diameter (m)
% rho - Density of water (kg/m^3)
% gamma - Surface tension of water (N/m)
%
% OUTPUT
% Td - surface temperature of particle (K)

lambda = free_path( T, D, M);
R = 8.3144621;
%x=T-10;
f= @(x)(x - T - ((D*M*L)/(R*ka))* ( (p/T)- (water_pvap(x)*kelvin_ratio( dp,x,rho,gamma,M )./x) )* FuchsSutugin( dp, lambda ));
Td= fzero(f,T);

end

