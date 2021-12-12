function V = value_function(k1,beta,alpha,delta,zgrid,piz,lss,kgrid,V0,i1,i2)

global c l

if k1 >= (1-delta)*kgrid(i1)
    
    EV = piz(i2,:)*transpose(interp1(kgrid,V0,k1,'linear','extrap'));
    
    labor_supply = @(l) (1-alpha)*zgrid(i2)*kgrid(i1)^(alpha)*l^(-1-alpha)+k1-zgrid(i2)*kgrid(i1)^(alpha)*l^(1-alpha)-(1-delta)*kgrid(i1);
    
    l = bisection(labor_supply,0,5*lss);
    
    c = zgrid(i2)*kgrid(i1)^(alpha)*l^(1-alpha)+(1-delta)*kgrid(i1)-k1;
    
    V = log(c)-0.5*l^2+beta*EV;
            
else
    
    V = -1e7;
    
end

V = -V;

end
    
    