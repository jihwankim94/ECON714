% ECON 714. Quant Macro-Econ Theory
% Problem Set 1


% 2. Integration

clear all
clc

format long g

A = 0;
B = 100;

rho = 0.04;
lambda = 0.02;

f = @(x) exp(-rho*x).*(-exp(-(1-exp(-lambda*x))));

nVals = [10 100 1000 10000 100000 1000000 10000000 100000000];
nrows = length(nVals);
ncols = 5;

solutions = zeros(nrows,ncols);
times = zeros(nrows,ncols);

for i = 1:nrows
   
    N = nVals(i);
    solutions(i,1) = N;
    times(i,1) = N;
    
    % Midpoint
    
    midpoint = 2;
    
    tic;
    
    for n = 1:N
        
        a = A+(n-1)*(B-A)/N;
        b = A+n*(B-A)/N;
    
        solutions(i,midpoint) = solutions(i,midpoint)+(b-a)*f((a+b)/2);
        
    end
            
    times(i,midpoint) = toc;
    
    % Trapezoid
    
    trapezoid = 3;
    
    tic;
    
    for n = 1:N
        
        a = A+(n-1)*(B-A)/N;
        b = A+n*(B-A)/N;
    
        solutions(i,trapezoid) = solutions(i,trapezoid)+(b-a)*(f(a)+f(b))/2;
        
    end
    
    times(i,trapezoid) = toc;
    
    % Simpson
    
    simpson = 4;
    
    tic;
    
    for n = 1:N
        
        a = A+(n-1)*(B-A)/N;
        b = A+n*(B-A)/N;
    
        solutions(i,simpson) = solutions(i,simpson)+((b-a)/6)*(f(a)+4*f((a+b)/2)+f(b));
        
    end
    
    times(i,simpson) = toc;
    
    % Monte Carlo
    
    montecarlo = 5;
    
    tic;
    
    rng(1)
    u = (B-A)*rand(N,1)+A;
    
    solutions(i,montecarlo) = (B-A)*mean(f(u));
    
    times(i,montecarlo) = toc;
        
end

%writematrix(solutions,'solutions.xls')
%writematrix(times,'times.xls')
