function [ dpdt ] = p_growth_rate( T, rho, Ntot, dp, ddpdt, M )
% Calculates the rate of change in vapour pressure during the growth
%
% INPUT
% T - Surrounding temperature (K)
% rho - Density of water (kg/m^3)
% Ntot - Number density of particles (#/m^3)
% dp - particle diameter (m)
% ddpdt - Rate of change in particle diameter (m/s)
% M - Molar mass of water (kg/mol)
%
% OUTPUT
% dpdt - Rate of change in partial pressure of water vapour (Pa/s)

R = 8.3144621;
dpdt = -(R*T*rho*Ntot*pi*dp^2*ddpdt)/(2*M);
end

