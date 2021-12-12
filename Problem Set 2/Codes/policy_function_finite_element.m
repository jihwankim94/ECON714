function [gc,gl,gk] = policy_function_finite_element(beta,alpha,delta,zgrid,Nz,piz,lss,Nk,kmin,kmax,kgrid,Nelement,psik,theta)

gc = zeros(Nk,Nz);
gl = zeros(Nk,Nz);
gk = zeros(Nk,Nz);

for i2 = 1:Nz
    
    for i1 = 1:Nk
    
    theta_l = theta((i2-1)*Nelement+1:i2*Nelement);
    
    gl(i1,i2) = dot(theta_l,psik(:,i1));
    
    y = zgrid(i2)*kgrid(i1)^(alpha)*gl(i1,i2)^(1-alpha);
    gc(i1,i2) = (1-alpha)*zgrid(i2)*kgrid(i1)^(alpha)*gl(i1,i2)^(-1-alpha);
    gk(i1,i2) = y+(1-delta)*kgrid(i1)-gc(i1,i2);
    
    end
    
end

end