function zero = price(p,x,w,sume)
% equilibrium prices
    
global n m

w_inverse = w.^(-1);

for i = 1:m
            
        c = zeros(1,n);
        d = zeros(n,1);
        
        for k = 1:n
            
            if i == 1
                
                c(k) = 1;
                d(k) = x(k,i)^(w(k,i)*w_inverse(k,i));
            
            else
                
                c(k) = p(i-1)^(w_inverse(k,i));
                d(k) = x(k,1)^(w(k,1)*w_inverse(k,i));
                
            end
            
        end
        
        zero(i) = sume(i)-c*d;
    
end

end