function [up1] = solveWaveEqn(u,up1,um1,t,dims,c,sloc,source)
   
    %Solving u_tt=c^2*uxx +f
    %L is the length(x), D the depth (y) and T the total time
    %c is the velocity, comes from background model BG
    %sloc  is the Source location
    %S  is the ricket wavelet (source)
    %I and V are the initial conditions.
    
    dt=dims.dt; 
    dx=dims.dx; %[m]
    dy=dims.dy; 
    Nt=dims.nt;
    Nx=dims.nx;
    Ny=dims.ny;
    dt2=dt^2;
    
    %Indices for slicing
    i=3:Nx-2; %indice in x direction
    j=3:Ny-2; %indice in y direction
 
    
    r2=(c.^2*dt^2)/dx^2; %Courant constant 
    beta=r2;
    alpha=(2-5*r2);
    
      
    
    % Inserting source    
    u(sloc)=u(sloc) + source(t,:);
    
    %solving the wave equation with a central differentiator
    %3th order: 9 point spacial stencil, 3 point time stencil
    
   
    
    up1(j,i) =  (beta(j,i)/12).*(-u(j-2,i)+16*u(j-1,i)-30*u(j,i)+16*u(j+1,i)-u(j+2,i) +...
                 -u(j,i-2)+16*u(j,i-1)-30*u(j,i)+16*u(j,i+1)-u(j,i+2)) +...
                 2*u(j,i)-um1(j,i);
end 
               
       
         
    
       
   
    