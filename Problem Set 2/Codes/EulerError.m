function e = EulerError(beta,alpha,delta,gc,gl,gk,zgrid,Nk,kgrid,piz,i2)

e = zeros(Nk,1);

for i1 = 1:Nk
    
    k = gk(i1,i2);
    c = interp1(kgrid,gc,k,'linear','extrap');
    l = interp1(kgrid,gl,k,'linear','extrap');
    
    e(i1,1) = log10(abs(1-gc(i1,i2)*beta*piz(i2,:)*(alpha*transpose(zgrid).*k^(alpha-1).*transpose(l).^(1-alpha).*transpose(c).^(-1))...
        -gc(i1,i2)*beta*piz(i2,:)*transpose(c).^(-1)*(1-delta))); 
    
        
end

end