function theta = theta_deep_learning(beta,alpha,delta,zgrid,Nz,piz,lss,Nk,kmin,kmax,kgrid,M,N,theta0)

options = optimset('Display','Iter','TolFun',1e-7,'TolX',1e-7);
theta = fsolve(@residual_finite_element,theta0,options);

function residual = residual_finite_element(theta0)

error = zeros(2*Nk,1);
residual = zeros(2*Nk*Nz,1);

theta0_V = theta0(1:N,1);
theta0_l = theta0(N+1:2*N,1);

rhs = zeros(Nz,1);

for i2 = 1:Nz
    
    V0 = zeros(Nk,1);
    gl = zeros(Nk,1);
    gk = zeros(Nk,1);
    gc = zeros(Nk,1);
    
    theta0_V_p = theta0_V((i2-1)*(5*M+1)+1:i2*(5*M+1));
    theta0_l_p = theta0_l((i2-1)*(5*M+1)+1:i2*(5*M+1));
    
    theta0_V01 = theta0_V_p(1:M,1);
    theta0_V11 = theta0_V_p(M+1:2*M,1);
    theta0_V02 = theta0_V_p(2*M+1:3*M,1);
    theta0_V12 = transpose(theta0_V_p(3*M+1:4*M,1));
    theta0_V03 = theta0_V_p(4*M+1,1);
    theta0_V13 = transpose(theta0_V_p(4*M+2:5*M+1,1));
 
    theta0_l01 = theta0_l_p(1:M,1);
    theta0_l11 = theta0_l_p(M+1:2*M,1);
    theta0_l02 = theta0_l_p(2*M+1:3*M,1);
    theta0_l12 = transpose(theta0_l_p(3*M+1:4*M,1));
    theta0_l03 = theta0_l_p(4*M+1,1);
    theta0_l13 = transpose(theta0_l_p(4*M+2:5*M+1,1));
    
    for i1 = 1:Nk
        
        k0 = kgrid(i1);
        
        zV1 = zeros(M,1);
        zV2 = zeros(M,1);
        phiV1 = zeros(M,1);
        phiV2 = zeros(M,1);
        
        zV1 = [theta0_V01 theta0_V11]*[1 ; k0];

%        phiV1 = (zV1>=0).*zV1;
        phiV1 = log(1+exp(zV1));
        
        for m = 1:M
            
        zV2(m) = theta0_V02(m)+theta0_V12*phiV1;
        
        end

%        phiV2 = (zV2>=0).*zV2;
        phiV2 = log(1+exp(zV2));

        V0(i1) = theta0_V03+theta0_V13*phiV2;
        
        zl1 = zeros(M,1);
        zl2 = zeros(M,1);
        phil1 = zeros(M,1);
        phil2 = zeros(M,1);
        
        zl1 = [theta0_l01 theta0_l11]*[1 ; k0];

%        phil1 = (zl1>=0).*zl1;
        phil1 = log(1+exp(zl1));
        
        for m = 1:M
            
            zl2(m) = theta0_l02(m)+theta0_l12*phil1;
        
        end

%        phil2 = (zl2>=0).*zl2;
        phil2 = log(1+exp(zl2));

        l0 = theta0_l03+theta0_l13*phil2;
        
        if l0 < 0
            
            l0 = 0;
            disp('l0 hit the lower bound')
            
        elseif l0 > 5*lss
            
            l0 = 5*lss;
            disp('l0 hit the upper bound')
            
        end
        
        gl(i1) = l0;
        
        y0 = zgrid(i2)*k0^(alpha)*l0^(1-alpha);
        c0 = max((1-alpha)*zgrid(i2)*k0^(alpha)*l0^(-1-alpha),1e-7);
        
        k1 = y0+(1-delta)*k0-c0;
        
        if k1 < 0.7*kmin
            
            k1 = 0.7*kmin;
            disp('k1 hit the lower bound')
            
        elseif k1 > 1.3*kmax
            
            k1 = 1.3*kmax;
            disp('k1 hit the upper bound')

        end
        
        gk(i1) = k1;
        gc(i1) = c0;
                
    end
        
    for i1 = 1:Nk
        
        V1 = zeros(Nz,1);
        
        for j2 = 1:Nz 
                
            theta0_V_p = theta0_V((j2-1)*(5*M+1)+1:j2*(5*M+1));
            theta0_l_p = theta0_l((j2-1)*(5*M+1)+1:j2*(5*M+1));
    
            theta0_V01 = theta0_V_p(1:M,1);
            theta0_V11 = theta0_V_p(M+1:2*M,1);
            theta0_V02 = theta0_V_p(2*M+1:3*M,1);
            theta0_V12 = transpose(theta0_V_p(3*M+1:4*M,1));
            theta0_V03 = theta0_V_p(4*M+1,1);
            theta0_V13 = transpose(theta0_V_p(4*M+2:5*M+1,1));
 
            theta0_l01 = theta0_l_p(1:M,1);
            theta0_l11 = theta0_l_p(M+1:2*M,1);
            theta0_l02 = theta0_l_p(2*M+1:3*M,1);
            theta0_l12 = transpose(theta0_l_p(3*M+1:4*M,1));
            theta0_l03 = theta0_l_p(4*M+1,1);
            theta0_l13 = transpose(theta0_l_p(4*M+2:5*M+1,1));
    
            zV1 = zeros(M,1);
            zV2 = zeros(M,1);
            phiV1 = zeros(M,1);
            phiV2 = zeros(M,1);
        
            zV1 = [theta0_V01 theta0_V11]*[1 ; k1];

%            phiV1 = (zV1>=0).*zV1;
            phiV1 = log(1+exp(zV1));

            for m = 1:M
            
                zV2(m) = theta0_V02(m)+theta0_V12*phiV1;
        
            end

%            phiV2 = (zV2>=0).*zV2;
            phiV2 = log(1+exp(zV2));

            V1(j2) = theta0_V03+theta0_V13*phiV2;
        
            zl1 = zeros(M,1);
            zl2 = zeros(M,1);
            phil1 = zeros(M,1);
            phil2 = zeros(M,1);
        
            zl1 = [theta0_l01 theta0_l11]*[1 ; k1];

%            phil1 = (zl1>=0).*zl1;
            phil1 = log(1+exp(zl1));
        
            for m = 1:M
            
                zl2(m) = theta0_l02(m)+theta0_l12*phil1;
        
            end

%            phil2 = (zl2>=0).*zl2;
            phil2 = log(1+exp(zl2));

            l1 = theta0_l03+theta0_l13*phil2;
            
            if l1 < 0
            
                l1 = 0;
                disp('l1 hit the lower bound')
            
            elseif l1 > 5*lss
            
                l1 = 5*lss;
                disp('l1 hit the upper bound')
            
            end
                
            y1 = zgrid(i2)*gk(i1)^(alpha)*l1^(1-alpha);
            c1 = max((1-alpha)*zgrid(i2)*gk(i1)^(alpha)*l1^(-1-alpha),1e-7);
        
            k2 = y1+(1-delta)*gk(i1)-c1;
        
            if k2 < 0.7*kmin
            
                k2 = 0.7*kmin;
                disp('k2 hit the lower bound')
            
            elseif k2 > 1.3*kmax
            
                k2 = 1.3*kmax;
                disp('k2 hit the upper bound')
                
            end
            
            rhs(j2,1) = 1/c1*(alpha*zgrid(j2)*gk(i1)^(alpha-1)*l1^(1-alpha)+1-delta);
            
        end
            
        EulerEquation_rhs = beta*dot(piz(i2,:),rhs);
        
        l0 = gl(i1);
        c0 = gc(i1);
        
        EulerEquation_lhs = 1/c0;
        
        BellmanEquation_rhs = log(c0)-0.5*l0^2+beta*dot(piz(i2,:),V1);
        BellmanEquation_lhs = V0(i1);
        
        error(i1) = EulerEquation_rhs-EulerEquation_lhs;
        error(Nk+i1) = BellmanEquation_rhs-BellmanEquation_lhs;
        
    end
        
    residual((i2-1)*Nk*2+1:i2*Nk*2) = error;

end
    
end

end
