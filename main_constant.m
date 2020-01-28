% MECH309 Final Project
% Main script

clear;
close all;

%% Constants and Parameters
M_inf = 0.85;
GAMMA = 1.4;
R = 287.058;    % J kg-1 K-1
T_inf = 293;    % K
P_inf = 100;    % kN/m^2
U_inf = M_inf*sqrt(GAMMA*R*T_inf);  % m/s

TOL = 1e-4;
MAXITER = 20000;
NUMPOINTS = 2001;

%% Calculations
% Airfoil
f_prime = @(x) 0.08*(-4*x+82);  % dy/dx

% Grid (coarse)
x = linspace(0,50,NUMPOINTS);
y = linspace(0,50,NUMPOINTS);
dx = x(2)-x(1);  
dy = y(2)-y(1);

% Initialize grid
phi = zeros(NUMPOINTS,NUMPOINTS);
c = 1/(dy^2);
f = 1/(dy^2);

% Initialize error
err = zeros(MAXITER,1);

tic
% Main loop
for n = 1:MAXITER           % each iteration
    % Store max error in each iteration
    maxErr = 0;
    % Update bottom boundary condition
    for i = 3:(0.4*(NUMPOINTS-1))                       % to the left of the airfoil
        phi(i,1) = phi(i,2);
    end
    for i = (0.4*(NUMPOINTS-1)+1):(0.42*(NUMPOINTS-1)+1)    % on the airfoil
        phi(i,1) = phi(i,2)-U_inf*feval(f_prime,x(i))*dy;
    end
    for i = (0.42*(NUMPOINTS-1)+1):(NUMPOINTS-1)        % to the right of the airfoil
        phi(i,1) = phi(i,2);
    end
    % Sweeping loop
    for i = 3:(NUMPOINTS-1)          % column by column
        for j = 2:(NUMPOINTS-1)
            phi_tmp = phi(i,j);
            % update Ai,j and Ai-i,j
            A1 = (1-M_inf^2)-(GAMMA+1)*M_inf^2/U_inf*(phi(i+1,j)-phi(i-1,j))/(2*dx);
            A2 = (1-M_inf^2)-(GAMMA+1)*M_inf^2/U_inf*(phi(i,j)-phi(i-2,j))/(2*dx);
            % update mui,j and mui-1,j according to A
            if A1>0
                mu1=0;
            else
                mu1=1;
            end
            if A2>0
                mu2=0;
            else
                mu2=1;
            end
            % update a, b, d, e
            a = (1-mu1)*A1/(dx^2);
            b = ((mu2*A2)-2*(1-mu1)*A1)/(dx^2)-2/(dy^2);
            d = ((1-mu1)*A1-2*mu2*A2)/(dx^2);
            e = (mu2*A2)/(dx^2);
            % Gauss-Seidal
            phi(i,j) = (-e*phi(i-2,j)-d*phi(i-1,j)-f*phi(i,j-1)...
                -c*phi(i,j+1)-a*phi(i+1,j))/b;
            % Calculate error
            maxErr = max(maxErr,abs(phi(i,j)-phi_tmp));
        end
    end
    % Store error
    err(n) = maxErr;
    disp(['Current iteration: ' num2str(n) ' | ' 'Current maxErr: ' num2str(maxErr)])
    % If error < tolerance, break
    if (err(n) < TOL)
        break;
    end
end
toc
loopTime = toc;

%% Convergence plot
figure(1);
semilogy(1:n,err(1:n));
title(['Convergence Plot, Mach Number=' num2str(M_inf)]);
xlabel('Number of iterations')
ylabel('L_{\infty}-norm error')

%% Pressure plot
cp = zeros(NUMPOINTS,NUMPOINTS);
for i = 2:NUMPOINTS-1
    for j = 1:NUMPOINTS
        cp(i,j) = -2/U_inf*(phi(i+1,j)-phi(i-1,j))/(2*dx);
    end
end
figure(2);
plot(x,-cp(:,1));
xlim([19.5 21.5])
title(['Coefficient of Pressure Plot, Mach Number=' num2str(M_inf)])
xlabel('x')
ylabel('C_{p}')

%% Pressure contour
u = zeros(NUMPOINTS,NUMPOINTS);
P = zeros(NUMPOINTS,NUMPOINTS);
for i = 2:NUMPOINTS-1
    for j = 1:NUMPOINTS
        u(i,j) = U_inf+(phi(i+1,j)-phi(i-1,j))/(2*dx);
        P(i,j) = P_inf*(1+(GAMMA-1)/2*M_inf^2*(1-(u(i,j)/U_inf)^2))^(GAMMA/(GAMMA-1));
    end
end
figure(3);
contour(x,y,P',50);
ylim([0 1]);
xlim([19.5 21.5]);
title(['Pressure Contour, Mach Number=' num2str(M_inf)])
xlabel('x')
ylabel('y')
colorbar
