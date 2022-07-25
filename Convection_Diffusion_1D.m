% 1D Convection-Diffusion Equation - Tt + u*Tx = a*Txx
% This code solves a specified 1D Partial Differential Convection-Diffusion
% Equation by using the methods FTCS or UPWIND according to user input

clc;
clear;

%% Specified Problem Data
a=10^-4;
u=0.05;
Dx=0.01; % Space step
xmin=0; xmax=1;
Nx=xmax/Dx; % Number of space steps
x=xmin:Dx:xmax;

%% Computation of Problem Parameters
while true
    prompt='Enter FTCS or UPWIND: ';
    answer=input(prompt,'s');
    if strcmp(answer,'FTCS')

        % FTCS METHOD
        Dt=0.01; % From the stability conditions C^2<=2s<=1 we find that 0<=Dt<=0.5. We choose Dt=0.1
        C=u*Dt/Dx;
        s=a*Dt/Dx^2;

        %Coefficients
        aC_e=-C/2;
        aD_e=s;
        aC_p=0;
        aD_p=-2*s;
        aC_w=C/2;
        aD_w=s;
        
        break;
        
    elseif strcmp(answer,'UPWIND')
        % UPWIND METHOD 
        Dt=0.01; % From the stability conditions C+2s<=1 we find that Dt<=0.475. We choose Dt=0.1
        C=u*Dt/Dx;
        s=a*Dt/Dx^2;

        % Coefficients
        if u>=0
            aC_e=0;
            aD_e=s;
            aC_p=-C;
            aD_p=-2*s;
            aC_w=C;
            aD_w=s;
        elseif u<0
            aC_e=-C;
            aD_e=s;
            aC_p=C;
            aD_p=-2*s;
            aC_w=0;
            aD_w=s;
        end
       
        break; 
    end
    disp('Incorrect input. Try again.')
end

ae=aC_e+aD_e;
ap=1+aC_p+aD_p;
aw=aC_w+aD_w;


%% Initialisation of matrix T 
% for the time we consider 100000 steps as we do not know how many steps are needed
T=zeros(Nx+1,10^5);

tolerance=10^-5;
error=1;
n=0; % iteration counter

%% Implementation of chosen method
while error>tolerance
    n=n+1;
    T(1,1)=10;  T(1,n+1)=10;  % Boundary condition T(0,t)=10
    T(Nx+1,1)=0; T(Nx+1,n+1)=0; % Boundary condition T(1,t)=0
     for j=2:Nx
         T(j,1)=0; % Initial condition T(x,0)=0

         T(j,n+1)=ae*T(j+1,n)+ap*T(j,n)+aw*T(j-1,n);
     end
    error=max(abs((T(2:Nx,n+1)-T(2:Nx,n))./T(2:Nx,n+1)));
end

fprintf('\n%s method \n',answer)
disp(['Number of iterations: ',num2str(n)])

%% Plotting
figure
plot(x,T(:,n+1)) % Final diagram
title(['T-x at ',num2str(n),' iterations'])
xlabel('x')
ylabel('T')

    