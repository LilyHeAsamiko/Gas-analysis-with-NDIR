function [lambda] = free_path( T, D, M)
% [lambda] = free_path( D, M, T ) Calculates free path length for molecule
% at temperature T(K) having diffusion coefficient D (m^2/s) and
% molar mass M (kg/mol)
%
% INPUT
% T - Surrounding temperature (K)
% D - Diffusion coefficient of water, (m^2/s)
% M - Molar mass of water (kg/mol)
%
% OUTPUT
% lambda - free path length (m)

R = 8.3144621;
lambda = 3.*D.*sqrt((pi.*M)./(8.*R.*T));


end

