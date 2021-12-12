function [V0,gc,gl,gk] = multigrid(beta,alpha,delta,zgrid,Nz,piz,lss,vss,kmin,kmax)

grid = [100 1000 10000];
Ngrid = length(grid);

for step = 1:Ngrid
   
    Nk = grid(step); % number of asset grids 
    
    kgrid = linspace(kmin,kmax,Nk); % asset grid
    
    V0 = ones(Nk,Nz)*vss; % Value function
    V1 = zeros(Nk,Nz);
    gc = zeros(Nk,Nz); % Consumption policy function
    gl = zeros(Nk,Nz); % Labor policy function
    gk = zeros(Nk,Nz); % Asset policy function
    
    if step > 1
        
        for i1 = 1:Nk
            
            V0(i1,:) = interp1(transpose(kgrid_previousstep),V_previousstep,kgrid(i1),'linear','extrap');
            
        end
        
    end
   
    distV = 1; 

    tol = 1e-7; % tolerance
    iter = 0; % number of iteration
    
    while distV > tol
    
       iter = iter+1;
       
       EV = piz*transpose(V0);
       
       for i2 = 1:Nz
           
           ik = 1;
           
           for i1 = 1:Nk
               
               V_updated = -1e7;
               
               for j1 = ik:Nk
                   
                   labor_supply = @(l) (1-alpha)*zgrid(i2)*kgrid(i1)^(alpha)*l^(-1-alpha)+kgrid(j1)-zgrid(i2)*kgrid(i1)^(alpha)*l^(1-alpha)-(1-delta)*kgrid(i1);
    
                   l = bisection(labor_supply,0,5*lss);
                   
                   c = zgrid(i2)*kgrid(i1)^(alpha)*l^(1-alpha)+(1-delta)*kgrid(i1)-kgrid(j1);
                   
                   V = log(c)-0.5*l^2+beta*EV(i2,j1);
                     
                   if V > V_updated
                       
                       V_updated = V;
                       consumption = c;
                       labor = l;
                       saving = kgrid(j1);
                       ik = j1;
                       
                   else 
                       
                       break;
                       
                   end
                   
               end
               
               V1(i1,i2) = V_updated;
               gc(i1,i2) = consumption;
               gl(i1,i2) = labor;
               gk(i1,i2) = saving;
               
           end
       
       end
       
       distV = max(max(abs(V1-V0)));
       V0 = V1;
       
    end
    
    if step < Ngrid
        
        V_previousstep = V0;
        kgrid_previousstep = kgrid;
    
    end
    
end

end
