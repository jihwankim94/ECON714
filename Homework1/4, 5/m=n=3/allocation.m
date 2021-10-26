function zero = allocation(x,alpha,w,lambda,sume)

global n m

for i = 1:m
    
    for j = 1:n
               
        w_inverse = w.^(-1);
       
        w_inverse(j,:) = [];
        
        alphalambda_ratio = (alpha.*lambda)*transpose((alpha.*lambda).^(-1));

        alphalambda_ratio(:,j) = [];
 
        a = zeros(1,n-1);
        b = zeros(n-1,1);
        
        for k = 1:n-1
            
            a(k) = alphalambda_ratio(j,k)^(w_inverse(k,i));
            b(k) = x(j,i)^(w(j,i)*w_inverse(k,i));
            
        end
        
        zero(j,i) = sume(i)-x(j,i)-a*b;

    end
    
end

end


