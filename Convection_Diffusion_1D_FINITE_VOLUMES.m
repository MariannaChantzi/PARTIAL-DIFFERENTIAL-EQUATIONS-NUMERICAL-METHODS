% 1D Convection-Diffusion Equation - Tt + u*Tx = a*Txx
% This code solves a specified 1D Partial Differential Convection-Diffusion Equation by
% using the FINITE VOLUMES method.

clc;
clear;

%% Specified Problem Data
k=10^-4;
u=-0.001;
xmin=0; xmax=1;
imax=21; %number of cell centers

%% Initialisation of total coefficients
aE=zeros(1,imax);
aP=zeros(1,imax);
aW=zeros(1,imax);

%% Grid generation
P=0.1;
Q=2;
[X,x,Dx,dxe,dxw]=mesh(imax,P,Q);

%% Coefficients for Diffusion
[aDE,aDP,aDW]=diffusion_coeff(k,dxe,dxw,imax);

%% Coefficients for Convection     
%UPWIND METHOD 
Dt=10^-2; %From the stability conditions C+2s<=1 we find that Dt<=0.02 for min(Dx). We choose Dt=0.01
[aCE,aCP,aCW]=convection_coeff_Upwind(u,imax);

%% Total coefficients
con=Dt./Dx; %constant matrix

aE=con.*(aDE-aCE);
aP=1+con.*(aDP-aCP);
aW=con.*(aDW-aCW);

%% Initialisation of matrix T 
%(for the time we consider 200000 steps as we do not know how many steps are needed).
T=zeros(1,imax);

%% Boundary conditions
T(1)=10;  %Boundary condition T(0,t)=0
T(imax)=0;  %Boundary condition T(1,t)=10

T_new=T;

%% Computation of T
tolerance=10^-5;
error=1;
n=0; %iteration counter

%% Time stamps for plots
n1=100/Dt;
n2=200/Dt;
n3=500/Dt;

while error>tolerance
    n=n+1;
     for j=2:imax-1
         T_new(j)=aE(j)*T(j+1)+aP(j)*T(j)+aW(j)*T(j-1);
     end
    error=max(abs((T_new(2:imax-1)-T(2:imax-1))./T_new(2:imax-1)));
    T=T_new;
    if n==n1
        figure(1)
        hold on
        plot(X,T_new(:),'b') %for t=100
    end
    if n==n2
         plot(X,T_new(:),'m') %for t=200
    end
    if n==n3
         plot(X,T_new(:),'g') %for t=500
    end
    
end

%% Plotting
xlabel('x')
ylabel('T')
title(sprintf('T-x using UPWIND method'))
if (n>n2)&&(n<n3)
    legend('t=100','t=200')
elseif (n>n1)&& (n<n2)
    legend('t=100')
elseif n<n1
else
    legend('t=100','t=200','t=500')
end
hold off

figure(2)
plot(X,T_new(:)) %final diagram
title(['F-x at ',num2str(n),' iterations'])
xlabel('x')
ylabel('T')

%% FUNCTIONS

function [X,x,Dx,dxe,dxw]=mesh(imax,P,Q)


%Initialisation
X=zeros(1,imax); %cell centers
x=zeros(1,imax); %cell faces
Dx=zeros(1,imax); %cell lengths
dxe=zeros(1,imax); %distance between cell centers in the east
dxw=zeros(1,imax); %distance between cell centers in the east

h=zeros(1,imax); %normalized variable
s=zeros(1,imax); %stretching function

xmin=0; xmax=1; %length of 1D grid is 1

dx=xmax/(imax-1); % imax=X(imax)/dx+1


for i=1:imax
    h(i)=((i-1)*dx-0)/(1-0);
    s(i)=P*h(i)+(1-P)*(1-tanh(Q*(1-h(i)))/tanh(Q));
    X(imax-(i-1))=xmax-s(i)*(xmax-xmin);
end

for i=1:imax-1
    x(i)=X(i)+(X(i+1)-X(i))/2;
end


for i=2:imax-1  
    Dx(i)=x(i)-x(i-1);
    dxe(i)=X(i+1)-X(i);
    dxw(i)=X(i)-X(i-1);
end
end

function [aDE,aDP,aDW]=diffusion_coeff(k,dxe,dxw,imax)

%Initialisation
    aDE=zeros(1,imax); aDP=zeros(1,imax); aDW=zeros(1,imax);
    
for i=2:imax-1
    aDE(i)=k/dxe(i);
    aDP(i)=-k/dxe(i)-k/dxw(i);
    aDW(i)=k/dxw(i);
end

end

function [aCE,aCP,aCW]=convection_coeff_Upwind(u,imax)

%Initialisation
    aCE=zeros(1,imax); aCP=zeros(1,imax); aCW=zeros(1,imax);
    
if u>=0
    for i=2:imax-1
        aCE(i)=0;
        aCP(i)=u;
        aCW(i)=-u;
    end
elseif u<0
    for i=2:imax-1
        aCE(i)=u;
        aCP(i)=-u;
        aCW(i)=0;
    end
end

end