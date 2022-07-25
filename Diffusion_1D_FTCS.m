% 1D Parabolic Diffusion Equation - Tt = a*Txx
% This code solves a specified 1D Partial Differential Diffusion Equation by using
% the explicit FTCS method. 

clc
clear

%% Specified Problem Data
a=10^(-4);
L=1; % Length of 1D direction
Dx=0.01; % Space step
Nx=L/Dx; % Number of space steps
Dt=0.1;  % Time step
Nt=20000; % Number of time steps 
s=a*Dt/(Dx^2); % s parameter

disp('The s parameter is')
disp(s)

%% Initialisation
T=zeros(Nx+1,Nt+1);
x=zeros(Nx+1,1);
t=zeros(Nt+1,1);

% The initial condition
for i=1:(Nx+1)
    T(i,1)=0;      % T(x,0)=0
    x(i)=(i-1)*Dx; % calculation of x vector
end

% The boundary conditions
for n=1:(Nt+1)
    T(1,n)=10;     % T(0,t)=10 
    T(Nx+1,n)=10;  % T(1,t)=10
    t(n)=(n-1)*Dt; % calculation of t vector
end

%% Implementation of the FTCS method
for n=1:Nt
    for i=2:Nx
        T(i,n+1)=s*T(i+1,n)+(1-2*s)*T(i,n)+s*T(i-1,n);
    end
end

%% Graphical representation at different selected times
figure(1)
plot(x,T(:,Nt/200+1))
title('Temperature at t1=10 (N1=100)')
xlabel('x')
ylabel('T')
hold on

plot(x,T(:,Nt/40+1))
title('Temperature at t2=50 (N2=500)')
xlabel('x')
ylabel('T')

plot(x,T(:,Nt/4+1))
title('Temperature at t3=500 (N3=5000)')
xlabel('x')
ylabel('T')
ylim([0 inf])

plot(x,T(:,Nt+1))
title('Temperature at t4=2000 (N4=20000)')
xlabel('x')
ylabel('T')
ylim([0 inf])
hold off








