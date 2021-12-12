% ECON 714. Quant Macro-Econ Theory
% Problem Set 2

clear all
clc

format longG

%%

% utility function

beta = 0.97;

% production function

alpha = 0.33;

delta = 0.1;

zgrid = [-0.05 0 0.05];
zgrid = exp(zgrid);
Nz = length(zgrid);

piz = [0.97 0.03 0.00 ; 0.01 0.98 0.01 ; 0.00 0.03 0.97];

%% 1. Steady State

l = @(k) ((1/beta-1+delta)/alpha)^(1/(1-alpha))*k;
c = @(k) k^(alpha)*l(k)^(1-alpha)-delta*k;

ss = @(k) 1/c(k)*(1-alpha)*k^(alpha)*l(k)^(-alpha)-l(k);

k0 = 1;
kss = fsolve(ss,k0)
lss = l(kss)
css = c(kss)

wss = (1-alpha)*kss^(alpha)*lss^(-alpha)
rss = alpha*kss^(alpha-1)*lss^(1-alpha)
yss = kss^(alpha)*lss^(1-alpha)

vss = log(css)-0.5*lss^2

%% 2. Value Function Iteration with a Fixed Grid

Nk = 250;

kmin = 0.7*kss;
kmax = 1.3*kss;

kgrid = linspace(kmin,kmax,Nk);

tic;

[V,gc,gl,gk] = fixedgrid(beta,alpha,delta,zgrid,Nz,piz,lss,vss,Nk,kgrid);

time = toc;

disp(' ')
disp(['Total Elapsed Time: ', num2str(time), ' seconds'])
disp(' ')

writematrix(gc,'gc_fixedgrid.xls')
writematrix(gl,'gl_fixedgrid.xls')
writematrix(gk,'gk_fixedgrid.xls')
writematrix(time,'time_fixedgrid.xls')

% Euler Error

EulerError_fixedgrid = zeros(Nk,Nz);

for i2 = 1:Nz
   
    EulerError_fixedgrid(:,i2) = EulerError(beta,alpha,delta,gc,gl,gk,zgrid,Nk,kgrid,piz,i2);
    
end

max_EulerError = max(max(EulerError_fixedgrid));
average_EulerError = sum(sum(EulerError_fixedgrid))/(Nk*Nz);

disp(' ')
disp(['Max Euler Error: ', num2str(max_EulerError), ', Average Euler Error: ', num2str(average_EulerError)])
disp(' ')

writematrix(max_EulerError,'maxEE_fixedgrid.xls')
writematrix(average_EulerError,'averageEE_fixedgrid.xls')

% plot

figure

plot(kgrid,gc(:,1),kgrid,gc(:,2),kgrid,gc(:,3))
title('Consumption Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Southeast')

figure

plot(kgrid,gl(:,1),kgrid,gl(:,2),kgrid,gl(:,3))
title('Labor Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Northeast')

figure

plot(kgrid,gk(:,1),kgrid,gk(:,2),kgrid,gk(:,3))
title('Capital Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Southeast')

%% 3. Accelerator

Nk = 250;

kmin = 0.7*kss;
kmax = 1.3*kss;

kgrid = linspace(kmin,kmax,Nk);

tic;

[V,gc,gl,gk] = accelerator(beta,alpha,delta,zgrid,Nz,piz,lss,vss,Nk,kgrid);

time = toc;

disp(' ')
disp(['Total Elapsed Time: ', num2str(time), ' seconds'])
disp(' ')

writematrix(gc,'gc_accelerator.xls')
writematrix(gl,'gl_accelerator.xls')
writematrix(gk,'gk_accelerator.xls')
writematrix(time,'time_accelerator.xls')

% Euler Error

EulerError_accelerator = zeros(Nk,Nz);

for i2 = 1:Nz
   
    EulerError_accelerator(:,i2) = EulerError(beta,alpha,delta,gc,gl,gk,zgrid,Nk,kgrid,piz,i2);
    
end

max_EulerError = max(max(EulerError_accelerator));
average_EulerError = sum(sum(EulerError_accelerator))/(Nk*Nz);

disp(' ')
disp(['Max Euler Error: ', num2str(max_EulerError), ', Average Euler Error: ', num2str(average_EulerError)])
disp(' ')

writematrix(max_EulerError,'maxEE_accelerator.xls')
writematrix(average_EulerError,'averageEE_accelerator.xls')

% plot

figure

plot(kgrid,gc(:,1),kgrid,gc(:,2),kgrid,gc(:,3))
title('Consumption Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Southeast')

figure

plot(kgrid,gl(:,1),kgrid,gl(:,2),kgrid,gl(:,3))
title('Labor Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Northeast')

figure

plot(kgrid,gk(:,1),kgrid,gk(:,2),kgrid,gk(:,3))
title('Capital Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Southeast')

%% 4. Multigrid

Nk = 10000;

kmin = 0.7*kss;
kmax = 1.3*kss;

kgrid = linspace(kmin,kmax,Nk);

tic;

[V,gc,gl,gk] = multigrid(beta,alpha,delta,zgrid,Nz,piz,lss,vss,kmin,kmax);

time = toc;

disp(' ')
disp(['Total Elapsed Time: ', num2str(time), ' seconds'])
disp(' ')

writematrix(gc,'gc_multigrid.xls')
writematrix(gl,'gl_multigrid.xls')
writematrix(gk,'gk_multigrid.xls')
writematrix(time,'time_multigrid.xls')

% Euler Error

EulerError_multigrid = zeros(Nk,Nz);

for i2 = 1:Nz
   
    EulerError_multigrid(:,i2) = EulerError(beta,alpha,delta,gc,gl,gk,zgrid,Nk,kgrid,piz,i2);
    
end

max_EulerError = max(max(EulerError_multigrid));
average_EulerError = sum(sum(EulerError_multigrid))/(Nk*Nz);

disp(' ')
disp(['Max Euler Error: ', num2str(max_EulerError), ', Average Euler Error: ', num2str(average_EulerError)])
disp(' ')

writematrix(max_EulerError,'maxEE_multigrid.xls')
writematrix(average_EulerError,'averageEE_multigrid.xls')

% plot

figure

plot(kgrid,V(:,1),kgrid,V(:,2),kgrid,V(:,3))
title('Value Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Southeast')

figure

plot(kgrid,gc(:,1),kgrid,gc(:,2),kgrid,gc(:,3))
title('Consumption Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Southeast')

figure

plot(kgrid,gl(:,1),kgrid,gl(:,2),kgrid,gl(:,3))
title('Labor Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Northeast')

figure

plot(kgrid,gk(:,1),kgrid,gk(:,2),kgrid,gk(:,3))
title('Capital Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Southeast')

%% 5. Chebyshev

Nk = 250;

kmin = 0.7*kss;
kmax = 1.3*kss;

kgrid = linspace(kmin,kmax,Nk);

Np = 5; % number of Chebyshev polynomials

N = Np*Nz;

tic;    
   
roots = -cos((2*transpose(1:Np)-1)*pi/(2*Np)); % roots of the Chebyshev polynomial

kgrid_chebyshev_1 = ((roots+1).*(kmax-kmin))./2+kmin;
    
T = ones(Np,Np); 
T(:,2) = roots;
    
for i2 = 3:Np
       
    T(:,i2) = 2.*roots.*T(:,i2-1)-T(:,i2-2);
        
end
    
theta0 = zeros(2*N,1);

for z = 1:Nz
            
    theta0((z-1)*Np+1) = vss;
    theta0((z-1)*Np+1+N) = lss;
            
end
        
theta = theta_chebyshev(beta,alpha,delta,zgrid,Nz,piz,lss,kmin,kmax,kgrid_chebyshev_1,Np,N,T,theta0);

[V,gc,gl,gk]= policy_function_chebyshev(beta,alpha,delta,Nz,zgrid,Nk,kmin,kmax,kgrid,Np,N,theta);

time = toc;

disp(' ')
disp(['Total Elapsed Time: ', num2str(time), ' seconds'])
disp(' ')

writematrix(gc,'gc_chebyshev.xls')
writematrix(gl,'gl_chebyshev.xls')
writematrix(gk,'gk_chebyshev.xls')
writematrix(time,'time_chebyshev.xls')

% Euler Error

EulerError_Chebyshev = zeros(Nk,Nz);

for i2 = 1:Nz
   
    EulerError_chebyshev(:,i2) = EulerError(beta,alpha,delta,gc,gl,gk,zgrid,Nk,kgrid,piz,i2);
    
end

max_EulerError = max(max(EulerError_chebyshev));
average_EulerError = sum(sum(EulerError_chebyshev))/(Nk*Nz);

disp(' ')
disp(['Max Euler Error: ', num2str(max_EulerError), ', Average Euler Error: ', num2str(average_EulerError)])
disp(' ')

writematrix(max_EulerError,'maxEE_chebyshev.xls')
writematrix(average_EulerError,'averageEE_chebyshev.xls')

% plot

figure

plot(kgrid,gc(:,1),kgrid,gc(:,2),kgrid,gc(:,3))
title('Consumption Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Southeast')

figure

plot(kgrid,gl(:,1),kgrid,gl(:,2),kgrid,gl(:,3))
title('Labor Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Northeast')

figure

plot(kgrid,gk(:,1),kgrid,gk(:,2),kgrid,gk(:,3))
title('Capital Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Southeast')

%% 6. Finite Elements

Nk = 250;

kmin = 0.7*kss;
kmax = 1.3*kss;

kgrid = linspace(kmin,kmax,Nk);

Nelement = 9;

kelement = linspace(kmin,kmax,Nelement);

tic

psi1 = @(x) (kelement(2)-x)/(kelement(2)-kelement(1)).*(kelement(1)<=x).*(x<kelement(2));
psi2 = @(x) (x-kelement(1))/(kelement(2)-kelement(1)).*(kelement(1)<=x).*(x<=kelement(2))+(kelement(3)-x)/(kelement(3)-kelement(2)).*(kelement(2)<x).*(x<kelement(3));
psi3 = @(x) (x-kelement(2))/(kelement(3)-kelement(2)).*(kelement(2)<=x).*(x<=kelement(3))+(kelement(4)-x)/(kelement(4)-kelement(3)).*(kelement(3)<x).*(x<kelement(4));
psi4 = @(x) (x-kelement(3))/(kelement(4)-kelement(3)).*(kelement(3)<=x).*(x<=kelement(4))+(kelement(5)-x)/(kelement(5)-kelement(4)).*(kelement(4)<x).*(x<kelement(5));
psi5 = @(x) (x-kelement(4))/(kelement(5)-kelement(4)).*(kelement(4)<=x).*(x<=kelement(5))+(kelement(6)-x)/(kelement(6)-kelement(5)).*(kelement(5)<x).*(x<kelement(6));
psi6 = @(x) (x-kelement(5))/(kelement(6)-kelement(5)).*(kelement(5)<=x).*(x<=kelement(6))+(kelement(7)-x)/(kelement(7)-kelement(6)).*(kelement(6)<x).*(x<kelement(7));
psi7 = @(x) (x-kelement(6))/(kelement(7)-kelement(6)).*(kelement(6)<=x).*(x<=kelement(7))+(kelement(8)-x)/(kelement(8)-kelement(7)).*(kelement(7)<x).*(x<kelement(8));
psi8 = @(x) (x-kelement(7))/(kelement(8)-kelement(7)).*(kelement(7)<=x).*(x<=kelement(8))+(kelement(9)-x)/(kelement(9)-kelement(8)).*(kelement(8)<x).*(x<kelement(9));
psi9 = @(x) (x-kelement(8))/(kelement(9)-kelement(8)).*(kelement(8)<=x).*(x<=kelement(9));

psi = @(x) [psi1(x) ; psi2(x) ; psi3(x) ; psi4(x) ; psi5(x) ; psi6(x) ; psi7(x) ; psi8(x) ; psi9(x)];

fplot(@(x) psi1(x),[kmin kmax])
hold on
fplot(@(x) psi2(x),[kmin kmax])
fplot(@(x) psi3(x),[kmin kmax])
fplot(@(x) psi4(x),[kmin kmax])
fplot(@(x) psi5(x),[kmin kmax])
fplot(@(x) psi6(x),[kmin kmax])
fplot(@(x) psi7(x),[kmin kmax])
fplot(@(x) psi8(x),[kmin kmax])
fplot(@(x) psi9(x),[kmin kmax])
hold off
title('Basis Functions')

psik1 = psi1(kgrid);
psik2 = psi2(kgrid);
psik3 = psi3(kgrid);
psik4 = psi4(kgrid);
psik5 = psi5(kgrid);
psik6 = psi6(kgrid);
psik7 = psi7(kgrid);
psik8 = psi8(kgrid);
psik9 = psi9(kgrid);

psik = [psik1 ; psik2 ; psik3 ; psik4 ; psik5 ; psik6 ; psik7 ; psik8 ; psik9];

psik_transpose = transpose(psik);

N = Nelement*Nz;

theta0 = lss*ones(Nelement*Nz,1);

theta = theta_finite_element(beta,alpha,delta,zgrid,Nz,piz,lss,Nk,kmin,kmax,kgrid,Nelement,psi,psik,psik_transpose,theta0)

[gc,gl,gk]= policy_function_finite_element(beta,alpha,delta,zgrid,Nz,piz,lss,Nk,kmin,kmax,kgrid,Nelement,psik,theta);

time = toc;

disp(' ')
disp(['Total Elapsed Time: ', num2str(time), ' seconds'])
disp(' ')

writematrix(gc,'gc_finite_element.xls')
writematrix(gl,'gl_finite_element.xls')
writematrix(gk,'gk_finite_element.xls')
writematrix(time,'time_finite_element.xls')

% Euler Error

EulerError_Chebyshev = zeros(Nk,Nz);

for i2 = 1:Nz
   
    EulerError_chebyshev(:,i2) = EulerError(beta,alpha,delta,gc,gl,gk,zgrid,Nk,kgrid,piz,i2);
    
end

max_EulerError = max(max(EulerError_chebyshev));
average_EulerError = sum(sum(EulerError_chebyshev))/(Nk*Nz);

disp(' ')
disp(['Max Euler Error: ', num2str(max_EulerError), ', Average Euler Error: ', num2str(average_EulerError)])
disp(' ')

writematrix(max_EulerError,'maxEE_finite_element.xls')
writematrix(average_EulerError,'averageEE_finite_element.xls')

% plot

figure

plot(kgrid,gc(:,1),kgrid,gc(:,2),kgrid,gc(:,3))
title('Consumption Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Southeast')

figure

plot(kgrid,gl(:,1),kgrid,gl(:,2),kgrid,gl(:,3))
title('Labor Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Northeast')

figure

plot(kgrid,gk(:,1),kgrid,gk(:,2),kgrid,gk(:,3))
title('Capital Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Southeast')

%% 7. Deep Learning

Nk = 250;

kmin = 0.7*kss;
kmax = 1.3*kss;

kgrid = linspace(kmin,kmax,Nk);

tic

M = 20; % number of nodes per layer

N = (5*M+1)*Nz;

theta0_V = ones(N,1)/(5*M+1);
theta0_V(M+1:2*M,1) = vss/kss;
theta0_V(6*M+2:7*M+1,1) = vss/kss;
theta0_V(11*M+3:12*M+2,1) = vss/kss;

theta0_l = ones(N,1)/(5*M+1);
theta0_l(M+1:2*M,1) = lss/kss;
theta0_l(6*M+2:7*M+1,1) = lss/kss;
theta0_l(11*M+3:12*M+2,1) = lss/kss;

%theta0_V = zeros(N,1);
%theta0_V(4*M+1,1) = vss;
%theta0_V(9*M+2,1) = vss;
%theta0_V(14*M+3,1) = vss;

%theta0_l = zeros(N,1);
%theta0_l(4*M+1,1) = lss;
%theta0_l(9*M+2,1) = lss;
%theta0_l(14*M+3,1) = lss;

theta0 = [theta0_V ; theta0_l];

theta = theta_deep_learning(beta,alpha,delta,zgrid,Nz,piz,lss,Nk,kmin,kmax,kgrid,M,N,theta0)

%%

[V,gc,gl,gk]= policy_function_deep_learning(beta,alpha,delta,zgrid,Nz,piz,lss,Nk,kmin,kmax,kgrid,M,N,theta);

time = toc;

disp(' ')
disp(['Total Elapsed Time: ', num2str(time), ' seconds'])
disp(' ')

%writematrix(gc,'gc_deep_learning.xls')
%writematrix(gl,'gl_deep_learning.xls')
%writematrix(gk,'gk_deep_learning.xls')
%writematrix(time,'time_deep_learning.xls')

% Euler Error

EulerError_Chebyshev = zeros(Nk,Nz);

for i2 = 1:Nz
   
    EulerError_chebyshev(:,i2) = EulerError(beta,alpha,delta,gc,gl,gk,zgrid,Nk,kgrid,piz,i2);
    
end

max_EulerError = max(max(EulerError_chebyshev));
average_EulerError = sum(sum(EulerError_chebyshev))/(Nk*Nz);

disp(' ')
disp(['Max Euler Error: ', num2str(max_EulerError), ', Average Euler Error: ', num2str(average_EulerError)])
disp(' ')

%writematrix(max_Eulererror,'maxEE_deep_learning.xls')
%writematrix(average_EulerError,'averageEE_deep_learning.xls')

% plot

figure

plot(kgrid,gc(:,1),kgrid,gc(:,2),kgrid,gc(:,3))
title('Consumption Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Southeast')

figure

plot(kgrid,gl(:,1),kgrid,gl(:,2),kgrid,gl(:,3))
title('Labor Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Northeast')

figure

plot(kgrid,gk(:,1),kgrid,gk(:,2),kgrid,gk(:,3))
title('Capital Function')
xlabel('Capital')
legend('z=-0.05','z=0','z=0.05','Location','Southeast')
