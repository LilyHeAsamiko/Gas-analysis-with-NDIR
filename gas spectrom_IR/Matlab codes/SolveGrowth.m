function [ t,dp,pw ] = SolveGrowth( T,D,M,L,ka,rho,gamma,Ntot,tmax,dp0,p0)
% This function solves particle condensation growth rate in case where 
% all the water vapour goes to the particles.
%
% INPUT
% T - Temperature (K)
% D - Diffusion coefficient of water, (m^2/s)
% M - Molar mass of water (kg/mol)
% L - Heat of evaporation of water (J/kg)
% ka - Thermal conductivity of air (W/(m*K))
% rho - Density of water (kg/m^3)
% gamma - Surface tension of water (N/m)
% Ntot - Number density of particles (#/m^3)
% tmax - How fas the calculation continues (s)
% dp0 - particle diameter in the beginning (m)
% p0 - partial pressure of water in the beginning (Pa)
%
% OUTPUT
% t - time(s)
% dp - particle diameter at the time t (m)
% pw - partial pressure of water at the time t (Pa)
%

options = odeset('RelTol',1e-6,'AbsTol',[1e-11 0.01],'InitialStep',1e-7,'Refine',10);
[t,Y] = ode15s(@(t,y) dyGrowth(t,y,T,D,M,L,ka,rho,gamma,Ntot,dp0),[0 tmax],[dp0 p0],options);
dp = Y(:,1);
pw = Y(:,2);
% figure
% 'stats','on'
% stream3(1:length(dp),1:length(pw),1:length(t),[dp,0,0],[0,pw,0],[0,0,t],dp0,p0,0);
end

