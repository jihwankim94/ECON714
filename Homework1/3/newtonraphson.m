function [X f] = newtonraphson(X0)
% Newton-Raphson

diffX = 1;
tol = 1e-5;

iter = 0;
maxiter = 1e5;

while iter < maxiter

    iter = iter+1;
    
    [f df hf] = rosenbrock(X0);

    X1 = X0-hf\df;
    
    diffX = abs(max(X1-X0));
            
    if diffX < tol
        
        break
        
    else
            
        X0 = X1;

    end
    
end

X = X1;

f = rosenbrock(X);

end
