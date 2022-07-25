% 1D Hyperbolic Convection Equation - Tt + u*Tx = 0
% This code solves a specified 1D Partial Differential Convection Equation by
% using the methods UPWIND, LAX-FRIEDRICHS, LAX-WENDROFF & FTCS

clc
clear

%% Specified Problem Data
u=1; xmax=1; 
Dx=0.01; % Space step
Nx=1/Dx; % Number of space steps
Dt=0.008; % Time step
Nt=100; % Number of time steps
C=u*Dt/Dx; % Courant constant

%% Initialisation
T=zeros(Nx+1,Nt+1);
x=linspace(0,1,Nx+1); %x vector

% Initial condition T(x,0)=0
for j=1:Nx+1
 T(j,1)=0;
end

% Boundary condition T(0,t)=10
for n=1:Nt+1
 T(1,n)=10;
end 

%% Upwind Method
for n=1:Nt
    for j=2:Nx+1
        T(j,n+1)=(1-C)*T(j,n)+C*T(j-1,n);
    end
end

figure(1)
plot(x,T(:,1:20:Nt+1))
title('Upwind method')
xlabel('x')
ylabel('T')
hold on

%% Lax-Friedrichs Method
for n=1:Nt
    for j=2:Nx
        T(j,n+1)=(1/2-C/2)*T(j+1,n)+(1/2+C/2)*T(j-1,n);
    end
end

figure(2)
plot(x,T(:,1:20:Nt+1))
title('Lax-Friedrichs method')
xlabel('x')
ylabel('T')

%% Lax-Wendroff Method
for n=1:Nt
    for j=2:Nx
        T(j,n+1)=(-C/2+(C^2)/2)*T(j+1,n)+(1-C^2)*T(j,n)+(C/2+(C^2)/2)*T(j-1,n);
    end
end

figure(3)
plot(x,T(:,1:20:Nt+1))
title('Lax-Wendroff method')
xlabel('x')
ylabel('T')


%% FTCS method
for n=1:Nt
    for j=2:Nx
        T(j,n+1)=T(j,n)+C/2*(T(j+1,n)-T(j-1,n));
    end
end

figure(4)
plot(x,T(:,1:20:Nt+1))
title('FTCS method')
xlabel('x')
ylabel('T')
hold off





