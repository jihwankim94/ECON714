function [X f] = steepestdescent(X0)

tol = 1e-5;

iter = 0;
maxiter = 1e5;

while iter < maxiter
        
    iter = iter+1;

    [f0 df0] =rosenbrock(X0);

    d = -df0; % search direction

    alpha = linesearch(@rosenbrock,X0,d,1e-5,1/8); % line search

    X1 = X0+alpha*d;
        
    [f1 df1] =rosenbrock(X1);

    if norm(df1) < tol
        
        break
         
    else
  
        X0 = X1;
       
    end
        
end

X = X1;

f = rosenbrock(X);

end