function [ ddpdt ] = dp_growth_rate( dp, T, rho, gamma, M, D, p, Td)
% This function calculates rate of change in particle diameter during the
% condensation growth.

% INPUT
% dp - particle diameter (m)
% T - Surrounding temperature (K)
% D - Diffusion coefficient of water, (m^2/s)
% M - Molar mass of water (kg/mol)
% rho - Density of water (kg/m^3)
% gamma - Surface tension of water (N/m)
% p - Surrounding partial pressure of water (Pa)
% Td - surface temperature of particle (K)
%
% OUTPUT
% ddpdt - Rate of change in particle diameter (m/s)


%Calculate the vapour pressure at the surface of the particle
pd = water_pvap( Td )*kelvin_ratio( dp,Td,rho,gamma,M );
lambda = free_path( T, D, M);
R = 8.3144621;

ddpdt = ( (4*D*M) / (R*rho*dp) ) * ( (p/T)- (pd/Td) ) * FuchsSutugin( dp, lambda );

end

