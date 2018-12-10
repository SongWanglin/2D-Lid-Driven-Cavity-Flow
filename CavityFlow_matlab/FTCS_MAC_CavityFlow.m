function FTCS_MAC()
%%%%%%%%%% Setup %%%%%%%%%%
% grid
Nc        = 256;         % cell number
Nx          = Nc + 2;    % cell number including ghost points
Ap = 1;               % aspect ratio
Lx = 1;    dx    = Lx/(Nx-1);  
Ly = Lx*Ap;            
dy = dx;   
Ny = floor(Ly/dy)+1;  
x  = 0:dx:Lx;                 y  = 0:dx:Ly;
% time step
dt        = 0.001;                % timestep
Final_Time = 2;   % final time
dtout      = 0.01;
Nsteps    = floor(Final_Time/dt);
nout      = floor(dtout/dt);    % number of time steps for output
%parameters
Re         = 500;       % Reynolds number
nu         = 1/Re;
cfl        = 0.5;       % cfl criterion
% allocate variables
p = ones(Nx,Ny);            u = zeros(Nx+1,Ny);        v = zeros(Nx,Ny+1);
Q = zeros(Nx-2,Ny-2);       F = zeros(Nx+1,Ny);        G = zeros(Nx,Ny+1);
% velocity profiles
un = 1;     vn = 0;    us = 0;     vs = 0; 
ue = 0;     ve = 0;    uw = 0;     vw = 0;
% boundary condtions %  
  % u-vel
  u(1,:) = 2*uw - u(2,:); u(end,:)= 2*ue-u(end-1,:);   % Left and rightwall
  u(:,1) = us;            u(:,end) = un;  % Bottom and top wall    
  % v-vel
  v(:,1) = 2*vs - v(:,2); v(:,end)= 2*vn - v(:,end-1);   % Bottom and top wall
  v(1,:) = vw;            v(end,:) = ve;   % Left and right wall
%%
%%%%%%%%%% main time loop %%%%%%%%%%
for tstep = 1:Nsteps
% check largest time step %
Umax = max(max(abs(u))); Vmax = max(max(abs(v)));
dt1 = 1/(Umax+Vmax)^2/cfl^2/Re;
dt2 = (0.25)*(dx^2)*Re;
dta = [dt1,dt2]; dt_max = min(dta);
if (dt > dt_max) , break;  end 
%%%%%%%%%% compute f,g,q %%%%%%%%%%
    i = 2:Nx-2;
    j = 2:Ny-1;
       F(i+1,j) = u(i+1,j) + dt*((u(i+2,  j) - 2.*u(i+1,j) + u(i,  j)) * nu/ (dx^2)...
         + (u(i+1,  j-1) - 2.*u(i+1,j) + u(i+1,  j+1)) * nu / ( dy^2)...
         - ((0.5.*(u(i+1,j) + u(i+2,j)) ).^2 - (0.5.*(u(i,j) + u(i+1,j)) ).^2)/dx...
         - ((u(i+1,j)+u(i+1,j+1)).*(v(i+1,j+1)+v(i,j+1))...
         - (u(i+1,j)+u(i+1,j-1)).*(v(i+1,j)+v(i,j)))*0.25/dy );
     
   i = 2:Nx-1;
    j = 2:Ny-2; 
        G(i,j+1) = v(i,j+1) + dt*((v(i+1,j+1) - 2*v(i,j+1) + v(i-1,j+1)) * nu/ ( dx^2)...
         +  (v(i,j+2) - 2*v(i,j+1) + v(i,j))* nu / ( dy^2)...
         - ((0.5.*(v(i,j+1) + v(i,j+2)) ).^2 - (0.5.*(v(i,j) + v(i,j+1)) ).^2)./dy...
         - ((u(i+1,j+1)+u(i+1,j)).*(v(i,j+1)+v(i+1,j+1))...
         - (u(i,j+1)+u(i,j)).*(v(i,j+1)+v(i-1,j+1)))*0.25/dx );
    % Q
    i = 2:Nx-1;
    j = 2:Ny-1;
    Q(i-1,j-1) = 1/dt*((F(i+1,j) - F(i,j))/dx + (G(i,j+1) - G(i,j))/dy);
    % Poisson solve
    p = SOR(Nx,Ny,dx,dy,u,v,p,Q,Re);
%%%%%%%%%% output %%%%%%%%%%
    if (mod(tstep,nout) == 0)
        time = tstep*dt;
        % Quiver Plot 
        [uplot, vplot ] = vel_plot(Nx,Ny,u,v);
        figure(1)
        clf,
        quiver(x,y,uplot,vplot,1.2,'k-')
        axis equal; axis([0 Lx 0 Ly])
        title({['2D Cavity Flow with Re = ',num2str(Re)];
            ['time(\itt) = ',num2str(time)]})
        xlabel('Spatial Coordinate (x) \rightarrow')
        ylabel('Spatial Coordinate (y) \rightarrow')
        drawnow;       
    end
    
%%%%%%%%%% update velocity %%%%%%%%%%
%%%%% u velocity
  i = 2:Nx-2;  j = 2:Ny-1;   
  u(i+1,j)    = F(i+1,j) - (dt/dx)*(p(i+1,j) - p(i,j));
  u(1,:) = 2*uw - u(2,:); u(end,:)= 2*ue-u(end-1,:);   % Left and rightwall
  u(:,1) = us;            u(:,end) = un;  % Bottom and top wall
%%%%% v velocity
  i = 2:Nx-1;  j = 2:Ny-2;
  v(i,j+1)    = G(i,j+1) - (dt/dy)*(p(i,j+1) - p(i,j));   
  v(:,1) = 2*vs - v(:,2); v(:,end)= 2*vn - v(:,end-1);   % Bottom and top wall
  v(1,:) = vw;            v(end,:) = ve;   % Left and right wall
end
 
%% plot the streamline, pressure contour and velocity field at final time
        figure(2)
        [uplot, vplot ] = vel_plot(Nx,Ny,u,v); skip = 2; 
        quiver(x(1:skip:end),y(1:skip:end),uplot(1:skip:end,1:skip:end),...
            vplot(1:skip:end,1:skip:end),1.5,'k-')
        axis equal; axis([0 Lx 0 Ly])
        title({['Velocity field for Re = ',num2str(Re)]})
        xlabel('Spatial Coordinate (x) \rightarrow')
        ylabel('Spatial Coordinate (y) \rightarrow')
        figure(3)
        sx = 0:.02:2;    sy = 0:.02:2;
        fn = stream2(x,y,uplot,vplot,sx,sy);
        fn1 = stream2(x,y,-uplot,-vplot,sx,sy);
        streamline(fn);  hold on;  streamline(fn1)
        axis equal; axis([0 Lx 0 Ly])
        title({['Streamline for Re = ',num2str(Re)]})
        xlabel('Spatial Coordinate (x) \rightarrow')
        ylabel('Spatial Coordinate (y) \rightarrow')
        figure(4)
        contourf(x,y,p',60,'w-');
        colormap(jet);  colorbar
        axis equal;axis([0 Lx 0 Ly])
        title({['Pressure Contour for Re = ',num2str(Re)]})
        xlabel('Spatial Coordinate (x) \rightarrow')
        ylabel('Spatial Coordinate (y) \rightarrow')
        %[yr,ur,xr,vr] = data_literature();    %if one side driven
        %[yr,ur,xr,vr] = data_literature2();  %if two side driven
        %[ur,yr,xr,vr] =analytical_Re1()   %if one side driven Re = 1
         %[ur,yr,xr,vr] =analytical_Re1()  %if one side driven A = 4/3
        figure(5)
        mid = Nx/2+1; 
        U(1:Nx,1:Ny)=(1/2)*(u(1:Nx,1:Ny)+ u(2:Nx+1,1:Ny));
        U_center = U(mid,:);plot(U_center,y)
        hold on;           % plot(ur,yr,'x')
        legend('program outcome','location','best');
        title({['Center u for Re = ',num2str(Re)]})
        xlabel('y');ylabel('Centerline u') ; hold off
        figure(6)
        mid = floor(Ny/2)+1; 
        V(1:Nx,1:Ny)=(1/2)*(v(1:Nx,1:Ny) + v(1:Nx,2:Ny+1));
        V_center = V(:,mid); 
        plot(x,V_center)
        legend('program outcome','location','best');
        title({['Center v for Re = ',num2str(Re)]})
        xlabel('x');ylabel('Centerline v') ; hold off
%  can use xlswrite to save u,v,p results for convergence study %       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [uplot,vplot] = vel_plot(Nx,Ny,u,v)
   uplot(1:Nx,1:Ny)  = 0.5 * (u(1:Nx,1:Ny) + u(2:Nx+1,1:Ny));
   vplot(1:Nx,1:Ny)  = 0.5 * (v(1:Nx,1:Ny) + v(1:Nx,2:Ny+1));
   Len   = sqrt(uplot.^2 + vplot.^2 + eps);
   uplot = (uplot./Len)';   vplot = (vplot./Len)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     poisson solve    %
function Successive_Over_Relaxation = SOR(Nx,Ny,dx,dy,u,v,p,Q,Re) 
    % Poisson matrix
    rf    = 1.6;      epsi   = 0.001;      itmax = 5000;
    change = 2*epsi;
    it     = 1;
    while (change > epsi)     
      pold = p;
    % boundary condition
        p(2:end-1 ,1    ) = p(2:end-1 ,2  ) - 2.0*v(2:end-1,3 )/(Re*dy); % bottom
        p(2:end-1 ,end  ) = p(2:end-1 ,end-1 ) + 2.0*v(2:end-1 ,end-2 )/(Re*dy); % top
        p(1   ,2:end-1) = p(2   ,2:end-1) - 2.0*u(3 ,2:end-1 )/(Re*dx); % left
        p(end ,2:end-1) = p(end-1 ,2:end-1) + 2.0*u(end-2 ,2:end-1 )/(Re*dx); % right    
      % SOR
        for i=2:Nx-1
            for j=2:Ny-1
                p(i,j) = 0.25*(p(i-1,j)+p(i,j-1)+p(i+1,j)+p(i,j+1) - ...
                            Q(i-1,j-1)*dx^2);
                p(i,j) = pold(i,j) + rf*(p(i,j)-pold(i,j));
            end
        end
        pmax = max(abs(pold(:)));
        if (pmax == 0) , pmax = 1.0;  end
        change =  max(abs( (pold(:)- p(:))/pmax ));
        it = it + 1;
        if (it > itmax) , break;  end   
        Successive_Over_Relaxation = p;
    end
