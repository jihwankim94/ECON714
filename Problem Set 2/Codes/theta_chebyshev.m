function theta = theta_chebyshev(beta,alpha,delta,zgrid,Nz,piz,lss,kmin,kmax,kgrid_chebyshev_1,Np,N,T,theta0)

options = optimset('Display','Iter','TolFun',1e-7,'TolX',1e-7);
theta = fsolve(@residual_chebyshev,theta0,options);

function residual = residual_chebyshev(theta0)

error = zeros(2*Np,1);
residual = zeros(2*N,1);

theta0_V = theta0(1:N,1);
theta0_l = theta0(N+1:2*N,1);

rhs = zeros(Nz,1);

for i2 = 1:Nz
    
    V0 = zeros(Np,1);
    gl = zeros(Np,1);
    gk = zeros(Np,1);
    gc = zeros(Np,1);
    
    theta0_V_p = theta0_V((i2-1)*Np+1:i2*Np);
    theta0_l_p = theta0_l((i2-1)*Np+1:i2*Np);
    
    for i1 = 1:Np
        
        V0(i1) = dot(theta0_V_p,T(i1,:));
        l0 = dot(theta0_l_p,T(i1,:));
        k0 = kgrid_chebyshev_1(i1);
        
        if l0 < 0
            
            l0 = 0
            
        elseif l0 > 5*lss
            
            l0 = 5*lss
                        
        end
        
        gl(i1) = l0;
        
        y0 = zgrid(i2)*k0^(alpha)*l0^(1-alpha);
        c0 = max((1-alpha)*zgrid(i2)*k0^(alpha)*l0^(-1-alpha),1e-7);
        
        k1 = y0+(1-delta)*k0-c0;
        
        if k1 < 0.7*kmin
            
            k1 = 0.7*kmin
            
        elseif k1 > 1.3*kmax
            
            k1 = 1.3*kmax
            
        end
        
        gk(i1) = k1;
        gc(i1) = c0;
        
    end
    
    kgrid_chebyshev_2 = (2*gk-(kmin+kmax))/(kmax-kmin);
    
    T_k = ones(Np,Np);
    T_k(:,2) = kgrid_chebyshev_2;
    
    for j2 = 3:Np
       
        T_k(:,j2) = 2*kgrid_chebyshev_2.*T_k(:,j2-1)-T_k(:,j2-2);
        
    end
    
    for i1 = 1:Np
        
        V1 = zeros(Nz,1);
        
        for j2 = 1:Nz
            
            theta0_V_p = theta0_V((j2-1)*Np+1:j2*Np);
            theta0_l_p = theta0_l((j2-1)*Np+1:j2*Np);
            V1(j2) = dot(theta0_V_p,T_k(i1,:));
            l1 = dot(theta0_l_p,T_k(i1,:));
            
            if l1 < 0
                
                l1 = 0;
                
            elseif l1 > 5*lss
                
                l1 = 5*lss;
                
            end

            y1 = zgrid(j2)*gk(i1)^(alpha)*l1^(1-alpha);
            c1 = max((1-alpha)*zgrid(j2)*gk(i1)^(alpha)*l1^(1-alpha),1e-7);
            
            k2 = y1+(1-delta)*gk(i1)-c1;
            
            if k2 < 0.7*kmin
                
                k2 = 0.7*kmin
                
            elseif k2 > 1.3*kmax
                
                k2 = 1.3*kmax
                
            end
            
            rhs(j2,1) = 1/c1*(alpha*zgrid(j2)*gk(i1)^(alpha-1)*l1^(-1-alpha)+1-delta);
            
        end
        
        EulerEquation_rhs = beta*dot(piz(i2,:),rhs);
        
        l0 = gl(i1);
        c0 = gc(i1);
                    
        EulerEquation_lhs = 1/c0;
        
        BellmanEquation_rhs = log(c0)-0.5*l0^2+beta*dot(piz(i2,:),V1);       
        BellmanEquation_lhs = V0(i1);
                    
        error(i1) = EulerEquation_rhs-EulerEquation_lhs;
        error(Np+i1) = BellmanEquation_rhs-BellmanEquation_lhs;
        
    end
    
    residual((i2-1)*Np*2+1:i2*Np*2) = error;

end
        
end
        
end                    