gnjfunction theta = theta_finite_element(beta,alpha,delta,zgrid,Nz,piz,lss,Nk,kmin,kmax,kgrid,Nelement,psi,psik,psik_transpose,theta0)

options = optimset('Display','Iter','TolFun',1e-7,'TolX',1e-7);
theta = fsolve(@galerkin,theta0,options)

function integral = galerkin(theta0)

error = zeros(Nk,1);
residual = zeros(Nk*Nz,1);

rhs = zeros(Nz,1);

for i2 = 1:Nz
    
    V0 = zeros(Nk,1);
    gl = zeros(Nk,1);
    gk = zeros(Nk,1);
    gc = zeros(Nk,1);
    
    theta0_l = theta0((i2-1)*Nelement+1:i2*Nelement);
    
    for i1 = 1:Nk
                
        k0 = kgrid(i1);
        l0 = dot(theta0_l,psik(:,i1));
        
        if l0 < 0
            
            l0 = 0;
            
        elseif l0 > 5*lss
            
            l0 = 5*lss;
            
        end
        
        gl(i1) = l0;
        
        y0 = zgrid(i2)*k0^(alpha)*l0^(1-alpha);
        c0 = max((1-alpha)*zgrid(i2)*k0^(alpha)*l0^(-1-alpha),1e-7);
        
        k1 = y0+(1-delta)*k0-c0;
        
        if k1 < 0.7*kmin
            
            k1 = 0.7*kmin;
            
        elseif k1 > 1.3*kmax
            
            k1 = 1.3*kmax;
            
        end
        
        gk(i1) = k1;
        gc(i1) = c0;
        
    end
    
    for i1 = 1:Nk
                
        for j2 = 1:Nz
            
            theta0_l = theta0((j2-1)*Nelement+1:j2*Nelement);
            l1 = dot(theta0_l,psi(gk(i1)));
            
            if l1 < 0
                
                l1 = 0;
                
            elseif l1 > 5*lss
                
                l1 = 5*lss;
                
            end

            y1 = zgrid(j2)*gk(i1)^(alpha)*l1^(1-alpha);
            c1 = max((1-alpha)*zgrid(j2)*gk(i1)^(alpha)*l1^(-1-alpha),1e-7);
            
            k2 = y1+(1-delta)*gk(i1)-c1;
            
            if k2 < 0.7*kmin
                
                k2 = 0.7*kmin;
                
            elseif k2 > 1.3*kmax
                
                k2 = 1.3*kmax;
                
            end
            
            rhs(j2,1) = 1/c1*(alpha*zgrid(j2)*gk(i1)^(alpha-1)*l1^(1-alpha)+1-delta);
            
        end
        
        EulerEquation_rhs = beta*dot(piz(i2,:),rhs);
        
        l0 = gl(i1);
        c0 = gc(i1);
                    
        EulerEquation_lhs = 1/c0;
        
        error(i1) = EulerEquation_rhs-EulerEquation_lhs;
        
    end
    
    residual((i2-1)*Nk+1:i2*Nk) = error;

end

residual = reshape(residual,[Nk,Nz]);  

residual_weighted = zeros(Nk,Nz*Nelement);

for j2 = 1:Nz

    for i2 = 1:Nelement
        
        residual_weighted(:,(i2-1)*Nz+j2) = psik_transpose(:,i2).*residual(:,j2);
        
    end

end

integral = sum(residual_weighted);

end

end