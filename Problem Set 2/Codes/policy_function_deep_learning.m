function [V,gc,gl,gk] = policy_function_deep_learning(beta,alpha,delta,zgrid,Nz,piz,lss,Nk,kmin,kmax,kgrid,M,N,theta)

V = zeros(Nk,Nz);
gc = zeros(Nk,Nz);
gl = zeros(Nk,Nz);
gk = zeros(Nk,Nz);

theta_V = theta(1:N,1);
theta_l = theta(N+1:2*N,1);

for i2 = 1:Nz
    
    for i1 = 1:Nk
    
    theta_V_p = theta_V((i2-1)*(5*M+1)+1:i2*(5*M+1));
    theta_l_p = theta_l((i2-1)*(5*M+1)+1:i2*(5*M+1));
    
    theta_V01 = theta_V_p(1:M,1);
    theta_V11 = theta_V_p(M+1:2*M,1);
    theta_V02 = theta_V_p(2*M+1:3*M,1);
    theta_V12 = transpose(theta_V_p(3*M+1:4*M,1));
    theta_V03 = theta_V_p(4*M+1,1);
    theta_V13 = transpose(theta_V_p(4*M+2:5*M+1,1));
 
    theta_l01 = theta_l_p(1:M,1);
    theta_l11 = theta_l_p(M+1:2*M,1);
    theta_l02 = theta_l_p(2*M+1:3*M,1);
    theta_l12 = transpose(theta_l_p(3*M+1:4*M,1));
    theta_l03 = theta_l_p(4*M+1,1);
    theta_l13 = transpose(theta_l_p(4*M+2:5*M+1,1));
            
    zV1 = zeros(M,1);
    zV2 = zeros(M,1);
    phiV1 = zeros(M,1);
    phiV2 = zeros(M,1);
        
    zV1 = [theta_V01 theta_V11]*[1 ; kgrid(i1)];

    phiV1 = (zV1>=0).*zV1;
        
    for m = 1:M
            
        zV2(m) = theta_V02(m)+theta_V12*phiV1;
        
    end

    phiV2 = (zV2>=0).*zV2;

    V(i1,i2) = theta_V03+theta_V13*phiV2;
        
    zl1 = zeros(M,1);
    zl2 = zeros(M,1);
    phil1 = zeros(M,1);
    phil2 = zeros(M,1);
        
    zl1 = [theta_l01 theta_l11]*[1 ; kgrid(i1)];

    phil1 = (zl1>=0).*zl1;
        
    for m = 1:M
            
        zl2(m) = theta_l02(m)+theta_l12*phil1;
        
    end

    phil2 = (zl2>=0).*zl2;

    gl(i1,i2) = theta_l03+theta_l13*phil2;
            
    y = zgrid(i2)*kgrid(i1)^(alpha)*gl(i1,i2)^(1-alpha);
    gc(i1,i2) = (1-alpha)*zgrid(i2)*kgrid(i1)^(alpha)*gl(i1,i2)^(-1-alpha);
    gk(i1,i2) = y+(1-delta)*kgrid(i1)-gc(i1,i2);
    
    end
    
end

end