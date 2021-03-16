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
nx = 11;
ny = 11;
nt = 100;                   % used during testing, no longer req
every = floor(nt / 1);      % plot frequency
dx = 1.0/(nx-1);
dy = 1.0/(ny-1);
sigma = .009;
nu=.01;
Re = 1/nu;
dt = sigma*dx*dy/nu; % CFL condition
rho = 1;
x = linspace(0,1,nx);
y = linspace(0,1,ny);
[x, y] = meshgrid(x,y);  
u = zeros(ny,nx); 
v = zeros(ny,nx);
un=zeros(ny,nx);
vn=zeros(ny,nx);
p = zeros(nx,ny);                   
varepsilon = 5e-7;
n=0;
delta = 1; % to enter the main loop
q = zeros(nx,ny); % normalized vector field
lefthole = floor(3*(ny-1)/4):ny-2;
righthole = 2:floor((ny-1)/4);
% Chorin projection method modification, using single intermediate u*
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
                    - dt ...vs
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

    u = us; 
    v = vs;

    if mod(n,every)==0 
        q = (u.^2 + v.^2); q = sqrt(q);     % vector field length
        unormed = (u./q); vnormed = (v./q); % normalizing vector field
%         contourf(x,y,p.',10,'w-'); hold on; % pressure field
%         colorbar;                 
%         caxis([min(min(p)) max(max(p))]);   % color axis            
%         contour(x,y,q.',10,'k-');  hold on; % velocity countour
        quiver(x,y,us.',vs.',1,'k');    % velocity field
        axis([-0.1 1.1 -0.1 1.1]); hold off;
        fprintf('%d. pmax = %f umax = %f vmax = %f delta = %.6f\n',...
            n,max(max(abs(p))),max(max(abs(u))),max(max(abs(v))),delta);
        colorbar vert; 
        
        drawnow;
    end

    n = n + 1;
end

