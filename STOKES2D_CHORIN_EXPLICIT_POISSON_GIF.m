%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Stokes 2D Chorin Projection method                       %
%                                                                         %
% Solving system of equations                                             %
% u_t + u * u_x + v * u_y = -1/rho * p_x + u_xx + u_yy                    %
% v_t + u * v_x + v * v_y = -1/rho * p_y + v_xx + v_yy                    %
% u_x + v_y = 0                                                           %
%                                                                         %
%                       by Rauan Kelesbekov                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Square 31x31 domain: 
% %%%%%%%%%%%%%%%%%%%%%%%%
% _______________________%
% inlet__________________%
% _______________________%
% %______________________%
% %______________________%
% %_______________________
% %_________________outlet
% %_______________________
% %%%%%%%%%%%%%%%%%%%%%%%%
% 
% IC: P = 0; u = 0; v = 0;
% 
% BC: 
% Walls: P_n  = 0; u = 0; v = 0;
% inlet BC: P=1; u = 1; v = 0;
% outlet BC: P=0; u_x=0; u=0;
% 
% Results obtained numerically via Chorin projection method:
% (1) Implicit burgers solution using upwind difference scheme. 
%       Preditcion step.
% (2) Over-Relaxation method for Poisson - pressure equation. 
%       Staggered mesh used in f_i for stability. 
% (3) Correction of step (1) by applying pressure field 
%       from (2) using staggered mesh. 
% 
% Result could be seen in the video below
% https://www.youtube.com/watch?v=xrbOS0AIqsg
% Velocity field is normalized (all vector lenghts are equal to 1)
% Velocity uncoloured black contour added for 10 different 
%   speeds from min to max
% Coloured pressure  contour added for 10 different
%   pressures from min to max
%  
% The code below also produces a GIF file namded 
% 'stokes2d_Chorin_Explicit_Overrelaxation.gif'

% variables declarations
clear; close all; clc;
format long;
figure('Renderer', 'painters', 'Position', [1000 600 850 800]);%plotcoord
nx = 31;
ny = 31;
nt = 100;                       % used during testing, no longer req
every = floor(nt / 1 * 10000);  % plot frequency
dx = 1.0/(nx-1);
dy = 1.0/(ny-1);
sigma = .009;
nu=.01;
Re = 1/nu;
dt = sigma*dx*dy/nu;            % CFL condition
rho = 1;
x = linspace(0,1,nx);
y = linspace(0,1,ny);
[x, y] = meshgrid(x,y);  
u = zeros(ny,nx); 
v = zeros(ny,nx);
un=zeros(ny,nx);
vn=zeros(ny,nx);
p = zeros(nx,ny);                   
varepsilon = 1e-7;
delta = 1;          % to enter the main loop
q = zeros(nx,ny);   % normalized vector field
lefthole = floor(3*(ny-1)/4):ny-5;
righthole = 4:floor((ny-1)/4);
% Chorin projection method modification, using single intermediate u*
figure(1);
filename = 'stokes2d_Chorin_Explicit_Overrelaxation.gif';
ctr = 1;
plotctr = 0;

while delta > varepsilon
    uprev = u; % velocity convergence parameters
    vprev = v;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   solving Burger's eq'n for u* v*                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    u(1,lefthole)=1;                        %inlet velocity Direchlet
    u(nx,righthole)=u(nx-1,righthole);      %outlet velocity Neumann
    % explicitly solving Burger's equation
    un=u;
    vn=v;
    for i=2:(ny-1)
        for j=2:(nx-1)
            % upwind difference scheme for stability
            un(i,j)= u(i,j) ...
                    - dt ...
                    * (max(u(i,j),0) * (u(i,j) - u(i-1,j))/dx ...
                    + min(u(i,j),0) * (u(i+1,j) - u(i,j))/dx)...
                    - dt ...
                    * (max(v(i,j),0) * (u(i,j) - u(i,j-1))/dy ...
                    + min(v(i,j),0) * (u(i,j+1) - u(i,j))/dx) ... 
                    + (1/Re*dt/dx^2) * (u(i+1,j)-2*u(i,j)+u(i-1,j)) ...
                    + (1/Re*dt/dy^2) * (u(i,j+1)-2*u(i,j)+u(i,j-1));
            vn(i,j)= v(i,j) ... 
                    - dt ...
                    * (max(u(i,j),0) * (v(i,j) - v(i-1,j))/dx ...
                    + min(u(i,j),0) * (v(i+1,j) - v(i,j))/dx)...
                    - dt ...
                    * (max(v(i,j),0) * (v(i,j) - v(i,j-1))/dy ...
                    + min(v(i,j),0) * (v(i,j+1) - v(i,j))/dx) ...
                    + (1/Re*dt/dx^2) *( v(i+1,j)-2*v(i,j)+v(i-1,j)) ...
                    + (1/Re*dt/dy^2) * (v(i,j+1)-2*v(i,j)+v(i,j-1));
        end
    end
    vs = vn; 
    us = un;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       solving Poisson eq                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m = 1; % pressure convergence parameter
    %boundary conditions
    p(1,:)=p(2,:);  
    p(nx,:)=p(nx-1,:);
    p(:,ny)=p(:,ny-1);
    p(:,1)=p(:,2);
    p(1,lefthole)=1;            % inlet pressure
    p(nx,righthole)=0;          % outlet pressure
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
                f(i,j) = - (rho/dt) ... 
                    * ((us(i+1,j) - us(i-1,j))/2/dx...
                    + (vs(i,j+1) - vs(i,j-1))/2/dy);
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
    delta = max(max(abs(uprev-u)))+max(max(abs(vprev-v)));
    u = uplus1;
    v = vplus1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Fancy plots                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plotctr = plotctr+delta;
    if mod(ctr,every)==1 || plotctr > 1e-1
        q = (u.^2 + v.^2); q = sqrt(q);     % vector field length
        unormed = (u./q); vnormed = (v./q); % normalizing vector field
        contourf(x,y,p.',20,'w-'); hold on; % pressure field
        colorbar;                 
        caxis([min(min(p)) max(max(p))]);   % color axis            
        contour(x,y,q.',20,'k-');  hold on; % velocity countour
        quiver(x,y,unormed.',vnormed.',0.4,'k');    % velocity field
        axis([-0.1 1.1 -0.1 1.1]); hold on;
        fprintf('%d. pmax = %f umax = %f vmax = %f delta = %.6f\n',...
            ctr,max(max(abs(p))),max(max(abs(u))),max(max(abs(v))),delta);
        colorbar vert; 
        drawnow;
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if ctr == 1
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
              imwrite(imind,cm,filename,'gif','WriteMode','append');
            end
            plotctr = 0;
    end
    ctr = ctr + 1;
end
