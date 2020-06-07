% Extinction efficiency calculation
a=(2+6)/2; % a=4;
m=1.33+0.001i; % refractive index of water
lambda=635*1e-9; % wavelength of light (nm)
%dp = logspace(1,4,1000); 
dp = [5 10 20]*a*1e-9;% particle vector, in which sizes the 
L=1;
N=1e10;
% calculate Mie-parameters: 
Qext = zeros(size(dp)); % Set extinction vector

for k = 1:numel(dp)
    result = Mie(m, pi*dp(k)/lambda); % Mie-parametrers
    Qext(k) = result(4); % 4. value of the vector is the extinction 
% efficiency, other results are not needed here.
end
% extinction efficiencies are calculated (nm)
y = exp(-pi*N*dp.^2.*Qext/4*1);

% Plot extinction efficiencies
figure(1);
loglog(dp,Qext,'-k')
xlabel('d_p (nm)')
ylabel('Q_{ext}')

figure(2)
plot(dp*1e9, y, '-k')
xlabel('d_p (nm)')
ylabel('I/I_{0}')

figure(3)
plot(t, y, '-k')
xlabel('t (s)')
ylabel('I/I_{0}')