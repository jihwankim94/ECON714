function [X f] = bfgs(X0)
% BFGS

h = eye(2); % hessian matrix

tol = 1e-5;

iter = 0;
maxiter = 1e5;

while iter < maxiter
    
    iter = iter+1;

    [f0 df0] =rosenbrock(X0);

    d = -h\df0; % search direction
    
    alpha = linesearch(@rosenbrock,X0,d,1e-5,1/8); % line search
    
    X1 = X0+alpha*d;
    
    [f1 df1] = rosenbrock(X1);

    s = X1-X0;
    y = df1-df0;
    
    h = h+(y*transpose(y))/(transpose(y)*s)-(h*s*transpose(s)*transpose(h))/(transpose(s)*h*s);
    
    if norm(df1) < tol

        break

    else

        X0 = X1;
        
    end
    
end

X = X1;

f = rosenbrock(X);

end