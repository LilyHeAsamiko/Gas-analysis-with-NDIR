% project 2
clear all
close all
%% question 5 supersaturation calculations
a=(2+6)/2; % a=4
dp=[5 10 20 50 100]*a*1e-9; % particle size (m)
%solve the supersaturation and plot how the function looks like
s=supersaturations(dp);
% error analysis
CI=zeros(5,2);
for i=1:5
    SEM=std(s(i))/sqrt(length(s));
    ts=tinv([0.025 0.975],length(s));
    CI(i,:)=s(i)+ts*SEM;
    display('CI for s',num2str(CI));
end

%'CI for s' 1.1124      1.0547       1.027      1.0107      1.0053 

%% question 6 and question 7
Pi=[991:0.001:1090];% customed initial pressure
d_p=particlesize(Pi);%
figure;
plot(1:length(Pi),d_p);
xlabel('initial pressure(mbar)');
ylabel('particle size (m)');
title('particle size v.s. initial pressure:991-1090')
%boxplot
i=int64(length(d_p)/10);
figure;
boxplot([d_p(1:1+i)',d_p(2+i:2+2*i)',d_p(3+2*i:3+3*i)',d_p(4+3*i:4+4*i)',d_p(5+4*i:5+5*i)',d_p(6+5*i:6+6*i)',d_p(7+6*i:7+7*i)',d_p(8+7*i:8+8*i)',d_p(9+8*i:9+9*i)'],'Notch','on')

D = 0.282e-4; %Diffusion coefficient of water, (m^2/s)
M = 18.016e-3; %Molar mass of water (kg/mol)
L = 2260e3; %Heat of evaporation of water (J/kg)
a=(6+2)/2; % a=4
A=10.23;
b=1750;
c=38;
Ti=296.15 %Tamperature (K)
rho = 1000; %Density of water (kg/m^3)
ka = 0.0257; %Thermal conductivity of air (W/(m*K))
Ntot = 1e10; %Number density of particles (#/m^3)
gamma = 72.8e-3; %Surface tension of water (N/m)
Pf=990 % (mbar)
dp=[5*a, 10*a, 20*a, 50*a, 100*a]*1e-9;
for i=1:length(dp)
Td(i)=droplet_temp(Tf,D,M,L,ka,Pf,dp(i),rho,gamma);
end
%solve the dfferences at the point s(a=4,dp=20)
difference=Pf-(10.^(A-b./(Td(1)-c)))*Pf./10.^(A-b./(Td-c));
%pi-pf
%(s(3)*water_pvap(Tf));
%Plot how the difference looks like
figure;
plot(dp,difference);
xlabel('particle size(m))');
ylabel('partial presure difference (mbar)');
title('partial presure difference (mbar) v.s. particle size ');
% a supersaturated state forms at cooling effect. When air tamperature 
% T goes down the saturated pressure ps decreases. Thus Pf is smaller than 
% Pi.  
for i=1:5
    SEM=std(difference(i))/sqrt(length(difference));
    ts=tinv([0.025 0.975],length(difference));
    CI(i)=difference(i)+ts*SEM;
end
display('CI for difference',num2str(CIU));
%'CI for difference   0    -67.25419     -172.4995      -357.252     -490.1837 =

 
%% question7
% solve the difference function numerically using interpolation
%for example, solve at the point dp=4*10^-10, 
dp=4*10^(-10);
P0=interp1(Pi,d_p,dp);%pf
%P0=1018.91100000000;
%error estimate
e=abs(particlesize(1018)-dp)/dp*100;
%e=0.0237;

%% question8
clear all;
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

tic
%Starting values for particle diameter (m) ja partial pressure of water (Pa). 
%Note, here saturation ratio in the beginning is 1.2.
dp0= 20*a*1e-9;
%p0 = 1.2*water_pvap(T);
d0=[5*a 10*a 20*a]*1e-9;
p0=supersaturations(d0)*water_pvap(T);
%Calculate and Plot the results
figure
subplot(1,2,1)
for i=1:3
[ t,dp,pw] = SolveGrowth( T,D,M,L,ka,rho,gamma,Ntot,tmax,dp0,p0(i) );
plot(t,dp)
hold on;
end
legend('5*a','10*a','20*a')
xlabel('t (s)')
ylabel('d_p (m)')

subplot(1,2,2)
for i=1:3
[ t,dp,pw ] = SolveGrowth( T,D,M,L,ka,rho,gamma,Ntot,tmax,dp0,p0(i) );
plot(t,pw)
hold on;
end
legend('5*a','10*a','20*a')
xlabel('t (s)')
ylabel('P_w (Pa)')
%in solve growth options = odeset('RelTol',1e-6,'AbsTol',[1e-11 0.01],'InitialStep',1e-7,'Refine',10,'stats','on');
%[t,Y] = ode15s(@(t,y) dyGrowth(t,y,T,D,M,L,ka,rho,gamma,Ntot,dp0),[0 tmax],[dp0 p0],options);
%For IVP ODE solvers,the choise of ReTol and Abtol are very important to the computational
%time  ,'RelTol',1e-6,'AbsTol',[1e-11 0.01]->2.590849s,
%'RelTol',1e-5,'AbsTol',[1e-10 0.01]->2.558769s
toc
%boxplot
i=int64(length(dp)/10);
figure;
boxplot([dp(1:1+i),dp(2+i:2+2*i),dp(3+2*i:3+3*i),dp(4+3*i:4+4*i),dp(5+4*i:5+5*i),dp(6+5*i:6+6*i),dp(7+6*i:7+7*i),dp(8+7*i:8+8*i),dp(9+8*i:9+9*i)],'Notch','on')
i=int64(length(pw)/10);
figure;
boxplot([pw(1:1+i),pw(2+i:2+2*i),pw(3+2*i:3+3*i),pw(4+3*i:4+4*i),pw(5+4*i:5+5*i),pw(6+5*i:6+6*i),pw(7+6*i:7+7*i),pw(8+7*i:8+8*i),pw(9+8*i:9+9*i)],'Notch','on')

%% question 9
dp=[5 10 20 50 100]*a*1e-9;
t=4;
[Qext,yext]=EXTINCTIONEFFICIENCY(dp,t);
% Plot extinction and extinction efficiencies
figure
loglog(dp*1e9,Qext,'-k');
xlabel('d_p (nm)')
ylabel('Q_{ext}')

figure
plot(dp*1e9, yext, '-k')
xlabel('d_p (nm)')
ylabel('I/I_{0}')

figure
plot(1:length(yext), yext, '-k')
xlabel('t (s)')
ylabel('I/I_{0}')

[x,y]=meshgrid(0:0.1:length(yext),0:0.1:length(Qext));
u=x;
v=y;
figure
quiver(x,y,u,v);
startx=yext(1):0.1:yext;
starty=Qext(1):0.1:Qext;
streamline(x,y,u,v,startx,starty);

%% question 10
%corelated to time
clear all;
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
dp0= 20*a*1e-9;
%p0 = 1.2*water_pvap(T);
d0=[5*a 10*a 20*a ]*1e-9;
p0=supersaturations(d0)*water_pvap(T);

for i=1:3
[ t,dp,pw] = SolveGrowth( T,D,M,L,ka,rho,gamma,Ntot,tmax,dp0,p0(i) );
[Qext,yext]=EXTINCTIONEFFICIENCY(dp,t);
figure
subplot (1,3,1)

loglog(dp*1e9,Qext,'-k')
xlabel('d_p (nm)')
ylabel('Q_{ext}')

subplot (1,3,2)
plot(dp*1e9, yext, '-k')
xlabel('d_p (nm)')
ylabel('I/I_{0}')

subplot (1,3,3)
plot(1:length(yext), yext, '-k')
xlabel('t (s)')
ylabel('I/I_{0}')

figure
[x,y]=meshgrid(0:0.5:length(yext),0:0.5:length(Qext));
u=x;
v=y;
figure
quiver(x,y,u,v);
startx=yext(1):1:yext;
starty=Qext(1):1:Qext;
streamline(x,y,u,v,startx,starty);
end

% 142 successful steps
% 4 failed attempts
% 244 function evaluations
% 2 partial derivatives
% 30 LU decompositions
% 238 solutions of linear systems
% 138 successful steps
% 1 failed attempts
% 264 function evaluations
% 2 partial derivatives
% 29 LU decompositions
% 258 solutions of linear systems
% 31 successful steps
% 3 failed attempts
% 58 function evaluations
% 3 partial derivatives
% 12 LU decompositions
% 47 solutions of linear systems
% 142 successful steps
% 4 failed attempts
% 244 function evaluations
% 2 partial derivatives
% 30 LU decompositions
% 238 solutions of linear systems
% 138 successful steps
% 1 failed attempts
% 264 function evaluations
% 2 partial derivatives
% 29 LU decompositions
% 258 solutions of linear systems
% 31 successful steps
% 3 failed attempts
% 58 function evaluations
% 3 partial derivatives
% 12 LU decompositions
% 47 solutions of linear systems


