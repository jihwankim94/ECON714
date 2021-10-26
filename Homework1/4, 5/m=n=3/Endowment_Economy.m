% ECON 714. Quant Macro-Econ Theory
% Problem Set 1


% 4. Computing Pareto efficient allocations

clear all
clc

format long g

global n m

m = 3; % goods
n = 3; % agents

% utility function

alpha = ones(n,1);
%alpha = rand(n,1);

w = -1*ones(n,m);
%w = -1*rand(n,m)-1;

% weight

lambda = 1/n*ones(n,1);
%lambda = rand(n,1);

% endowment

e = ones(n,m);
%e = rand(n,m);
sume = sum(e);


% computing pareto efficient allocations

myfun1 = @(x) allocation(x,alpha,w,lambda,sume)


x0 = zeros(n,m);

for j = 1:n
    
    for i = 1:m
       
        x0(j,i) = 1/n*sume(i);
        
    end
    
end

x1 = fsolve(myfun1,x0) 


% 5. Computing Equilibrium allocations

myfun2 = @(p) price(p,x1,w,sume);

p0 = ones(1,m-1);

p1 = fsolve(myfun2,p0); 

p1 = [1 p1]

