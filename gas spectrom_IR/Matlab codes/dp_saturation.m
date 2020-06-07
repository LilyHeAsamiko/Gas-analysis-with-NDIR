
T = 296.15; %Temperature (K)
D = 0.282e-4; %Diffusion coefficient of water, (m^2/s)
M = 18.016e-3; %Molar mass of water (kg/mol)
L = 2260e3; %Heat of evaporation of water (J/kg)
ka = 0.0257; %Thermal conductivity of air (W/(m*K))
rho = 1000; %Density of water (kg/m^3)
gamma = 72.8e-3; %Surface tension of water (N/m)
Ntot = 1e10; %Number density of particles (#/m^3)
tmax = 2; %How fast the calculation continues (s)
a=(2+6)/2; %4

%Starting values for particle diameter (m) ja partial pressure of water (Pa). 
%Note, here saturation ratio in the beginning is 1.2.
dp0 = 20*a*1e-9;
%p0 = 1.2*water_pvap(T);
p0=supersaturations(5*a*1e-9)*water_pvap(T);

%Solve the growth of particles 
[ t,dp,pw ] = SolveGrowth( T,D,M,L,ka,rho,gamma,Ntot,tmax,dp0,p0 );


%Plot the results. More info "doc subplot".
subplot(1,2,1)
plot(t,dp,'-k')
xlabel('t (s)')
ylabel('d_p (m)')

subplot(1,2,2)
plot(t,pw,'-k')
xlabel('t (s)')
ylabel('P_w (Pa)')