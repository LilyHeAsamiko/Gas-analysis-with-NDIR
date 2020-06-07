%particle size calculations
function d_p=particlesize(Pi)
%parameters for water 17-1 
Ti=296.15;
Pf=990; %(mbar) room pressure
gamma= 72.8e-3; %Surface tension of water (N/m)
M = 18.016e-3; %Molar mass of water (kg/mol)
Ntot = 1e10; %Number density of particles (#/m^3)
v=M*Ntot*1e-17;
a=10.23;
b=1750;
c=38;
Tf=Ti*((Pf./Pi).^(gamma-1)./gamma);%17-13 use the Pf/Pi to calculate the Ti
psi=log(a-b./(Ti-c));%17-1 calculated the initial pressure at the supersaturation
psf=log(a-b./(Tf-c));%17-1 calculated the final pressure at the supersaturation
SR=psi./psf*Pf./Pi;%17-16 calculate the SR
d_p=4*v*gamma./log(SR);%17-3 calculate the d_p
end
