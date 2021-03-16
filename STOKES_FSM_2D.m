%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Stokes 2D Chorin Projection method                       %
%                                                                         %
% Solving system of equations                                             %
% u_t + u * u_x + v * u_y = -1/rho * p_x + u_xx + u_yy                    %
% v_t + u * v_x + v * v_y = -1/rho * p_y + v_xx + v_yy                    %
% u_x + v_y = 0                                                           %
% by Rauan Kelesbekov                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% variable declarations
clear; close all; clc;
format long;
figure('Renderer', 'painters', 'Position', [1000 600 850 800]);
fprintf('initialization')
nx = 51;
ny = 51;
nt = 40;
Lx = 1;
Ly = 1;
dx = Lx/(nx-1);
dy = Ly/(ny-1);
sigma = .009;
nu=.01;
Re = 1/nu;
dt = sigma*dx*dy/nu;    % CFL condition
dt = 0.1;               % Using stable FSM scheme now. no need for CFL      
rho = 1;
Re = 50; 
x = linspace(0,1,nx);
y = linspace(0,1,ny);
[x, y] = meshgrid(x,y);  
u = zeros(ny,nx); 
v = zeros(ny,nx);
un=zeros(ny,nx);
vn=zeros(ny,nx);
uplus1 = u;
vplus1 = v;
p = zeros(nx,ny);                   
varepsilon = 5e-7;
tol = 1e-6;
delta = 1; % to enter the main loop
q = zeros(nx,ny); % normalized vector field
inlet = floor(3*(ny-1)/4):ny-5;
outlet = 4:floor((ny-1)/4);
% Chorin projection method modification, using single intermediate u*
figure(1);
filename = 'stokes2d_FSM.gif';
ctr = 1;
plotctr = 0;
% boundary conditions
u(2:nx-1,1)=0;
u(1,2:ny-1)=0;
u(1,2:ny-1)=0;
u(nx,2:ny-1)=0;
u(2:nx-1,ny)=0;
v(2:nx-1,1)=0;
v(1,2:ny-1)=0;
v(ny,2:ny-1)=0;
v(2:nx-1,ny)=0;
u(1,inlet)=1;                        %inlet velocity Direchlet
u(nx,outlet)=u(nx-1,outlet);      %outlet velocity Neumann
fprintf(', time loop\n--20%%--40%%--60%%--80%%--100%%\n')
for k = 1:nt% while delta > varepsilon
%     u(2:nx-1,1)=0;
%     u(1,2:ny-1)=0;
%     u(1,2:ny-1)=0;
%     u(nx,2:ny-1)=0;
%     u(2:nx-1,ny)=0;
%     v(2:nx-1,1)=0;
%     v(1,2:ny-1)=0;
%     v(ny,2:ny-1)=0;
%     v(2:nx-1,ny)=0;
%     u(1,inlet)=1;                        %inlet velocity Direchlet
%     u(nx,outlet)=u(nx-1,outlet);      %outlet velocity Neumann
    uprev = u; % velocity convergence parameters
    vprev = v;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   solving Burger's eq'n for u* v*                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[us, vs] = FSMonce(Lx,Ly,nx,ny,dt,Re,varepsilon,inlet,outlet,u,v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       solving Poisson eq                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    m = 1; % pressure convergence parameter
    %boundary conditions
    p(1,:)=p(2,:);  
    p(nx,:)=p(nx-1,:);
    p(:,ny)=p(:,ny-1);
    p(:,1)=p(:,2);
    pmax = 1;
    p(1,inlet)=pmax;            % inlet pressure
    p(nx,outlet)=0;          % outlet pressure
    pplus1  = p;
    % calculating free term 
    f = zeros(nx,ny);
    for j = 2:ny-1
        for i = 2:nx-1
            usplus05    = (us(i+1,j) + us(i,j))/2;
            usminus05   = (us(i,j) + us(i-1,j))/2;
            vsplus05    = (vs(i,j+1) + vs(i,j))/2;
            vsminus05   = (vs(i,j) + vs(i,j-1))/2;
            f(i,j) = - (rho/dt) ...
                * ((usplus05-usminus05)/dx ...
                    + (vsplus05-vsminus05)/dy);
%                 f(i,j) = - (rho/dt) ... 
%                     * ((us(i+1,j) - us(i-1,j))/2/dx...
%                     + (vs(i,j+1) - vs(i,j-1))/2/dy);
        end
    end
    % relaxation method, solving Poisson equation, finding pressure field
    w = 1.95;   % optimal relaxation parameter
    while m > varepsilon
        p = pplus1;
        for i = 2:nx-1
            for j = 2:ny-1
                pplus1(i,j) = (w/4)*(p(i+1,j) + p(i,j+1) ...
                    - 4*(1-1/w)*p(i,j) + dx^2 * f(i,j)...
                    + pplus1(i-1,j) + pplus1(i,j-1));
            end
        end
        m = max(max(abs(pplus1) - abs(p)));
    end
    p = pplus1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    correcting u,v considering pressure                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    uplus1 = us;
    vplus1 = vs;
    for i = 2:nx-1 
        for j = 2:ny-1 
            % staggered mesh
            pxplus005   = (p(i+1,j) + p(i,j))/2;
            pxminus05   = (p(i,j) + p(i-1,j))/2;                        
            pyplus005   = (p(i,j+1) + p(i,j))/2;
            pyminus05   = (p(i,j) + p(i,j-1))/2;
            uplus1(i,j) = - dt / rho * (pxplus005-pxminus05)/dx + us(i,j);
            vplus1(i,j) = - dt / rho * (pyplus005-pyminus05)/dy + vs(i,j);
        end
    end
    uprev = u; % velocity convergence parameters
    vprev = v;
    delta = max(max(abs(uprev(3:nx-2,2:ny-2)-u(3:nx-2,2:ny-2))))...
        +max(max(abs(vprev(3:nx-2,2:ny-2)-v(3:nx-2,2:ny-2))));
    u = uplus1;
    v = vplus1;
%     u(nx,outlet)=u(nx-1,outlet);                                            % NO NEED FOR THIS TERM, SHALL CALCULATE WITHOUT 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Fancy plots                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ctr = ctr + 1;
%     if floor(25*k/nt)>floor(25*(k-1)/nt) || k == 1
        fprintf('.');
        q = (u.^2 + v.^2); q = sqrt(q);     % vector field length
        unormed = (u./q); vnormed = (v./q); % normalizing vector field
%         unormed = u; vnormed = v;
        contourf(x,y,p.',10,'w-'); hold on; % pressure field
        colorbar;                 
        caxis([min(min(p)) max(max(p))]);   % color axis            
        contour(x,y,q.',10,'k-');  hold on; % velocity countour
        quiver(x,y,unormed.',vnormed.',0.8,'k');    % velocity field
%         quiver(x,y,u.',v.',0.8,'k');
        axis([-0.1 1.1 -0.1 1.1]); hold off;
        colorbar vert; 
        drawnow;
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if k == 1
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
              imwrite(imind,cm,filename,'gif','WriteMode','append');
            end
%     end

end
fprintf('\n');
q = (u.^2 + v.^2); q = sqrt(q);     % vector field length
unormed = (u./q); vnormed = (v./q); % normalizing vector field
%         unormed = u; vnormed = v;
contourf(x,y,p.',20,'w-'); hold on; % pressure field
colorbar;                 
caxis([min(min(p)) max(max(p))]);   % color axis            
contour(x,y,q.',20,'k-');  hold on; % velocity countour
quiver(x,y,unormed.',vnormed.',0.8,'k');    % velocity field
axis([-0.1 1.1 -0.1 1.1]); hold off;
colorbar vert; 
drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if k == 1
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
      imwrite(imind,cm,filename,'gif','WriteMode','append');
    end