% question 5
% supersaturation calculations 
function s=supersaturations(d)
% Starting parameters
rho = 1000; %Density of water (kg/m^3)
M = 18.016e-3; %Molar mass of water (kg/mol)
T=296.15; % Room Tamperature (K)
v=M/rho; % molar volume of the liquid
gamma=72.8e-3; % surface tension
R=8.31446; %gas constant(J/K mol)

s=exp(4*gamma.*v./R./T./d);%17-3 for insoluble particles on which pure water condensates.
figure;
plot(d,s,'o-');
xlabel('particle size (m)')
ylabel('required supersaturation')
end
