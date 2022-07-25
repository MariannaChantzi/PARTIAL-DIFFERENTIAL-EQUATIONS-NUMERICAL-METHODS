% 2D POISSON EQUATION - Txx + Tyy = S(x,y)
% This code solves a specified 2D Partial Differential Poisson Equation by
% using the methods Jacobi, Gauss-Seidel & S.O.R (for various relaxation coefficients).

clear;
clc;

%% Specified Problem Data & Initialisation
Dx=0.02; Dy=0.02; % Space step
x=0:Dx:1;         % Vector x
y=0:Dy:1;         % Vector y
jmax=length(x);
kmax=length(y);
S=-1;             % Value of source in Poisson equation
tolerance=1.e-5;  % Maximum allowed value of error

%% Coefficients 
% T(j,k)=L*S(j,k)-M*(T(j+1,k)+T(j-1,k))-N*(T(j,k+1)+T(j,k-1))
L=1/(-2/Dx^2-2/Dy^2);
M=L/Dx^2;
N=L/Dy^2;


%% Jacobi method
% Initialisation of matrix and boundary conditions
T_J=zeros(jmax,kmax); % T(0,y)=T(1,y)=T(x,0)=T(x,1)=0
Told=T_J;
error=1;

N_J=0; % iteration counter
while error>tolerance
    N_J=N_J+1;
    for j=2:(jmax-1)
        for k=2:(kmax-1)
            T_J(j,k)=L*S-M*(Told(j+1,k)+Told(j-1,k))-N*(Told(j,k+1)+Told(j,k-1));
        end
    end
    error=max(abs((T_J(2:(jmax-1),2:(kmax-1))-Told(2:(jmax-1),2:(kmax-1)))./T_J(2:(jmax-1),2:(kmax-1))));
    Told=T_J;
end

disp('Jacobi Method:')
disp(['Number of iterations: ',num2str(N_J)])
figure(1)
surf(x,y,T_J)
title('Jacobi Method')
xlabel('x')
ylabel('y')
zlabel('T')


%% Gauss-Seidel method
T_GS=zeros(jmax,kmax);
Told=T_GS;
error=1;

N_GS=0; % iteration counter
while error>tolerance
    N_GS=N_GS+1;
    for j=2:(jmax-1)
        for k=2:(kmax-1)
            T_GS(j,k)=L*S-M*(Told(j+1,k)+T_GS(j-1,k))-N*(Told(j,k+1)+T_GS(j,k-1));
        end
    end
    error=max(abs((T_GS(2:(jmax-1),2:(kmax-1))-Told(2:(jmax-1),2:(kmax-1)))./T_GS(2:(jmax-1),2:(kmax-1))));
    Told=T_GS;
end

disp('Gauss-Seidel Method:')
disp(['Number of iterations: ',num2str(N_GS)])
figure(2)
surf(x,y,T_GS)
title('Gauss-Seidel Method')
xlabel('x')
ylabel('y')
zlabel('T')


%% Successive Over-Relaxation (S.O.R.) method

%--------------------------------------------
% For ë=0.5
T_SOR_1=zeros(jmax,kmax);
Told=T_SOR_1;
error=1;
lamda1=0.5;

N_SOR_1=0;
while error>tolerance
    N_SOR_1=N_SOR_1+1;
    for j=2:jmax-1
        for k=2:kmax-1
            T_SOR_1(j,k)=(1-lamda1)*Told(j,k)+lamda1*(L*S-M*(Told(j+1,k)+...
                +T_SOR_1(j-1,k))-N*(Told(j,k+1)+T_SOR_1(j,k-1)));
        end
    end
    error=max(abs((T_SOR_1(2:(jmax-1),2:(kmax-1))-Told(2:(jmax-1),2:(kmax-1)))./T_SOR_1(2:(jmax-1),2:(kmax-1))));
    Told=T_SOR_1;     
end

disp('S.O.R. Method (ë=0.5):')
disp(['Number of iterations (ë=0.5): ',num2str(N_SOR_1)])
figure(3)
surf(x,y,T_SOR_1)
title('S.O.R. Method (ë=0.5)')
xlabel('x')
ylabel('y')
zlabel('T')

%-------------------------------------------
% For ë=1.3
T_SOR_2=zeros(jmax,kmax);
Told=T_SOR_2;
error=1;
lamda2=1.3;

N_SOR_2=0;
while error>tolerance
    N_SOR_2=N_SOR_2+1;
    for j=2:jmax-1
        for k=2:kmax-1
            T_SOR_2(j,k)=(1-lamda2)*Told(j,k)+lamda2*(L*S-M*(Told(j+1,k)+...
                +T_SOR_2(j-1,k))-N*(Told(j,k+1)+T_SOR_2(j,k-1)));
        end
    end
    error=max(abs((T_SOR_2(2:(jmax-1),2:(kmax-1))-Told(2:(jmax-1),2:(kmax-1)))./T_SOR_2(2:(jmax-1),2:(kmax-1))));
    Told=T_SOR_2;     
end

disp('S.O.R. Method (ë=1.3):')
disp(['Number of iterations (ë=1.3): ',num2str(N_SOR_2)])
figure(4)
surf(x,y,T_SOR_2)
title('S.O.R. Method (ë=1.3)')
xlabel('x')
ylabel('y')
zlabel('T')

%-------------------------------------------
% For ë=1.8
T_SOR_3=zeros(jmax,kmax);
Told=T_SOR_3;
error=1;
lamda3=1.8;

N_SOR_3=0;
while error>tolerance
    N_SOR_3=N_SOR_3+1;
    for j=2:jmax-1
        for k=2:kmax-1
            T_SOR_3(j,k)=(1-lamda3)*Told(j,k)+lamda3*(L*S-M*(Told(j+1,k)+...
                +T_SOR_3(j-1,k))-N*(Told(j,k+1)+T_SOR_3(j,k-1)));
        end
    end
    error=max(abs((T_SOR_3(2:(jmax-1),2:(kmax-1))-Told(2:(jmax-1),2:(kmax-1)))./T_SOR_3(2:(jmax-1),2:(kmax-1))));
    Told=T_SOR_3;     
end

disp('S.O.R. Method (ë=1.8):')
disp(['Number of iterations (ë=1.8): ',num2str(N_SOR_3)])
figure(5)
surf(x,y,T_SOR_3)
title('S.O.R. Method (ë=1.8)')
xlabel('x')
ylabel('y')
zlabel('T')

