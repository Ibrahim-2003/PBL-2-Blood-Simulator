clc;
clear all;
close all force;
%% We want to solve the system of three first-order differential equations
% For that, first of all, let's list some of the parameters summarized on
% his book, as an example. We choose "25x25 TW", for example
alpha=0.01;
FMExchangeStifness=1.04*(10^(-11));
K0Anisotropy=0.85*(10^5);
KAnisotropy= 0.075*(10^5);
gamma=2.21*(10^5);
mu0=4*pi*(10^(-7));
Ms=10000*(1000/(4*pi));
HK=2*KAnisotropy/(mu0*Ms);
AppliedField=60*(10^(-3))*(10000*(1000/(4*pi)));
aParameter=(pi^2)/12;
% Now, we can define the time array
time=linspace(0,200,10000).*(10^(-9));
% Now, we can define the initial conditions
x0=0; % dimensionless
phi0=pi/2; % rad
Delta0=10*(10^(-9)); % nm
initial_conditions=[x0 phi0 Delta0];
%% Now, let's write our system of coupled differential equations
% First of all, let's create our system of differential equations
% symbolically
syms alphasym gammasym Hasym HKsym mu0sym Mssym asym Asym K0sym Ksym t x(t) y(t) z(t) T Y 
Dx=diff(x);
Dy=diff(y);
Dz=diff(z);
Eqn1=alphasym*Dx/z+Dy==gammasym*Hasym;
Eqn2=Dx/z-alphasym*Dy==gammasym*HKsym*sin(2*y)/2; 
Eqn3=Dz==(gammasym/(alphasym*mu0sym*Mssym*asym))*(Asym/z-(K0sym+Ksym*((sin(2*y))^2))*z);
Eqn1s=simplify(lhs(Eqn1)-rhs(Eqn1),'Steps',100);
Eqn2s=simplify(lhs(Eqn2)-rhs(Eqn2),'Steps',100);
Eqn3s=simplify(lhs(Eqn3)-rhs(Eqn3),'Steps',100);
[VF,Subs]=odeToVectorField(Eqn1s,Eqn2s,Eqn3s);
odefcn=matlabFunction(VF,'Vars',{T Y alphasym gammasym Hasym HKsym mu0sym Mssym asym Asym K0sym Ksym});
% Now, let's solve numerically the system of differential equations
odefcn=@(T,Y,alphasym,gammasym,Hasym,HKsym,mu0sym,Mssym,asym,Asym,K0sym,Ksym)...
    [(Y(3)./(alphasym.^2+1.0)).*(alphasym.*gammasym.*Hasym+gammasym.*HKsym.*sin(Y(2).*2.0)./2.0);...
    gammasym.*(Hasym-(alphasym./(1.0+alphasym.^2)).*(alphasym.*Hasym+HKsym.*sin(Y(2).*2.0)./2.0));...
    (gammasym./(alphasym.*mu0sym.*Mssym.*asym)).*(Asym./Y(3)-(K0sym+Ksym.*((sin(Y(2))).^2)).*Y(3))];
[Time,Variables]=ode45(@(T,Y)odefcn(T,Y,alpha,gamma,AppliedField,HK,mu0,Ms,aParameter,FMExchangeStifness,...
    K0Anisotropy,KAnisotropy),time,initial_conditions);

figure
for k = 1:size(Variables,2)
    subplot(size(Variables,2), 1, k)
    plot(Time,Variables(:,k))
    grid
    title(sprintf('Variable %d',k))
end

