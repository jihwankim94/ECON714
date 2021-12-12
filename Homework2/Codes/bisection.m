function xm = bisection(f,xl,xu)

distf = 1;

tol = 1e-4; 

iter = 0;
maxiter = 1e3;

if iter == maxiter

    disp('The maximum number of iteration is reached')

end

if f(xl)*f(xu)>0
    
    disp('Error')
    xm = NaN;
    
else
    
    iter = iter+1;
    
    while (distf > tol || iter > maxiter)
        
        xm = 0.5*xl+0.5*xu;
        
        if f(xm)>0
            
            xl = xm;
            
        else
            
            xu = xm;
            
        end
        
        distf = abs(f(xm));
        
    end
    
end

end