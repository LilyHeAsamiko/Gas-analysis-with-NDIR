function [ kr ] = kelvin_ratio( dp,T,rho,gamma,M )
% [ kr ] = kelvin_ratio( dp,T,rho,gamma,M ) Calculates kelvin ratio
%
% INPUT
% dp - particle diameter (m)
% T - Surrounding temperature (K)
% M - Molar mass of water (kg/mol)
% rho - Density of water (kg/m^3)
% gamma - Surface tension of water (N/m)
%
% OUTPUT
% lambda - free path length (m)

R = 8.3144621;
kr = exp((4.*M.*gamma)./(rho.*R.*T.*dp));


end

