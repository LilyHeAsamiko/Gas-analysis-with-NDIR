% Extinction efficiency calculation
function [Qext,ext]=EXTINCTIONEFFICIENCY(dp,t)
m=1.33+0.001i; % refractive index of water
lambda=635*1e-9;% wavelength of light (m)
%dp = logspace(1,4,1000); % particle vector, in which sizes the 
%dp = [5 10 20]*a*1e9;    % extinction efficiencies are calculated (m)
L=1; % length of tube
N=1e10;
% calculate Mie-parameters:
Qext = zeros(size(dp)); % Set extinction vector

for k = 1:numel(dp)
    result = Mie(m, pi*dp(k)/lambda); % Mie-parametrers
    Qext(k) = result(4); % 4. value of the vector is the extinction 
            % efficiency, other results are not needed here.
end
% extinction is calculated (m)
ext = exp(-pi*N*dp.^2.*Qext/4*L);
end