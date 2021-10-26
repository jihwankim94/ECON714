% ECON 714. Quant Macro-Econ Theory
% Homework 1
% 3. Optimization: basic problem

% https://academics.uccs.edu/rcascava/Math4420/codes.html

clear all
clc

format long g

X0 = [10 ; 10];

solutions = zeros(3,4);
times = zeros(1,4);

% Newton-Raphson

tic;

[X f] = newtonraphson(X0)

solutions(1,1) = X(1);
solutions(2,1) = X(2);
solutions(3,1) = f;

times(1,1) = toc;

% BFGS

tic;

[X f] = bfgs(X0)

solutions(1,2) = X(1);
solutions(2,2) = X(2);
solutions(3,2) = f;

times(1,2) = toc;

% Steepest Descent

tic;

[X f] = steepestdescent(X0)

solutions(1,3) = X(1);
solutions(2,3) = X(2);
solutions(3,3) = f;

times(1,3) = toc;

% Conjugate Gradient

tic;

[X f] = conjugategradient(X0)

solutions(1,4) = X(1);
solutions(2,4) = X(2);
solutions(3,4) = f;

times(1,4) = toc;

%writematrix(solutions,'solutions.xls')
%writematrix(times,'times.xls')
