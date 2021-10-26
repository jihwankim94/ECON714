function [X f] = conjugategradient(X0)

[f0 df0] = rosenbrock(X0);
    
d0 = -df0; % search direction
r0 = -df0;

tol = 1e-5;

iter = 0;
maxiter = 1e5;

while iter < maxiter

    iter = iter+1;

    alpha = linesearch(@rosenbrock,X0,d0,1e-5,1/8); % line search

    X1 = X0+alpha*d0;

    [f1 df1] = rosenbrock(X1);

    r1 = -df1;
    
    beta = (transpose(r1)*r1)/(transpose(r0)*r0);

    d1 = r1+beta*d0;

    if norm(df1) < tol
        
        break
         
    else
  
        X0 = X1;
        d0 = d1;
        r0 = r1;
       
    end

end

X = X1;

f = rosenbrock(X);

end