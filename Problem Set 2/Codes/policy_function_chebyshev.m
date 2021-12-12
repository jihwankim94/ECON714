function [V,gc,gl,gk] = policy_function_chebyshev(beta,alpha,delta,Nz,zgrid,Nk,kmin,kmax,kgrid,Np,N,theta)

kgrid_chebyshev_3 = (2*transpose(kgrid)-(kmin+kmax))/(kmax-kmin);

T_k = ones(Nk,Np);
T_k(:,2) = kgrid_chebyshev_3;

for i2 = 3:Np
    
   T_k(:,i2) = 2*kgrid_chebyshev_3.*T_k(:,i2-1)-T_k(:,i2-2);
    
end

theta_V = theta(1:N,1);
theta_l = theta(N+1:2*N,1);

V = zeros(Nk,Nz);
gc = zeros(Nk,Nz);
gl = zeros(Nk,Nz);
gk = zeros(Nk,Nz);

for i2 = 1:Nz
    
    for i1 = 1:Nk
    
    theta_V_p = theta_V((i2-1)*Np+1:i2*Np);
    theta_l_p = theta_l((i2-1)*Np+1:i2*Np);
    
    V(i1,i2) = dot(theta_V_p,T_k(i1,:));
    gl(i1,i2) = dot(theta_l_p,T_k(i1,:));
    
    y = zgrid(i2)*kgrid(i1)^(alpha)*gl(i1,i2)^(1-alpha);
    gc(i1,i2) = (1-alpha)*zgrid(i2)*kgrid(i1)^(alpha)*gl(i1,i2)^(-1-alpha);
    gk(i1,i2) = y+(1-delta)*kgrid(i1)-gc(i1,i2);
    
    end
    
end

end