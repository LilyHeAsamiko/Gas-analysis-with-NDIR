%usthe particle size to calculate the supersaturation with interpolation
%create growing particle size vector
a=(6+2)/2; % a=4
A=10.23;
b=1750;
c=38;
Tf=296.15; % initial tamperature (K);
rho = 1000; %Density of water (kg/m^3)
gamma = 72.8e-3; %Surface tension of water (N/m)
dp=[5*a, 10*a, 20*a, 50*a, 100*a]*1e-9;
pf=990 %(mbar)
for i=1:length(dp)
Td(i)=droplet_temp(Tf,D,M,L,ka,pf,dp(i),rho,gamma);
end
%solve the supersaturation and plot how the function looks like
s=supersaturations(dp);
%solve the dfferences at the point s(a=4,dp=20)
difference=pf-(10^(A-b/(T-c)))*pf./(s(3)*water_pvap(Td));
%Plot how the function looks like
figure;
plot(dp,difference);
% In when T goes down


