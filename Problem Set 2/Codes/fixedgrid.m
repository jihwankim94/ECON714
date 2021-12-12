function [V0,gc,gl,gk] = fixedgrid(beta,alpha,delta,zgrid,Nz,piz,lss,vss,Nk,kgrid)

global c l

V0 = ones(Nk,Nz)*vss; % Value function
V1 = zeros(Nk,Nz);
gc = zeros(Nk,Nz); % Consumption policy function
gl = zeros(Nk,Nz); % Labor policy function
gk = zeros(Nk,Nz); % Asset policy function

distV = 1; 

tol = 1e-7; % tolerance
iter = 0; % number of iteration

while distV > tol
    
    iter = iter+1;
    
    for i2 = 1:Nz
        
        for i1 = 1:Nk
       
            [k,V] = fminbnd(@(k) value_function(k,beta,alpha,delta,zgrid,piz,lss,kgrid,V0,i1,i2),0.7*kgrid(i1),1.3*kgrid(i1));
            
            V1(i1,i2) = -V;
            gc(i1,i2) = c;
            gl(i1,i2) = l;
            gk(i1,i2) = k;
            
        end
        
    end
    
    distV = max(max(abs(V1-V0)))
    V0 = V1;
    
end
            
end