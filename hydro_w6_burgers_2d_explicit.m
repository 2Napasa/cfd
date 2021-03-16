% 2D burgers eq'n 
% given  system
% u_t + u*u_x + v*u_y= (1/Re) * (u_xx+u_yy)
% v_t + u*v_x + v*v_y= (1/Re) * (v_xx+v_yy)
% initial condition zero everywhere for u and v
% boundary left side u(3:ny-2,1)=1; v(3:ny-2,1)=1; i.e. having positive x
% and y direction of velocity on left edge
% since the edges could cause problems first and last 2 nodes are zerofied.


% the solution displays the infinitisimal thickness flat pipe with the stream 
%   going to the North East direction
% considered viscosity (via Reynolds number), no Pressure and thus whole
%   fluid stream has nonchanging direction
% the discretized first derivative had central difference scheme, similarly
%   to 1D Burgers equation where the step forward backward to be used
%   dependent on the direction of the velocity, central difference scheme
%   is said to be absolutely stable, regardless of the direction of the
%   flow
% The plot attached below after the code. 

clear; close all; clc;
nx = 101;
ny = 101;
dx = 2.0/(nx-1);
dy = 2.0/(ny-1);
% CFL condition
sigma = .009; 
nu=0.01;
dt = sigma*dx*dy/nu;
% domain  discretization
x = linspace(0,1,nx);
y = linspace(0,1,ny);
%i.c.
un=ones(ny,nx)*0;
vn=ones(ny,nx)*0;
 

[x, y] = meshgrid(x,y);                            
u=un+1;
v=vn+1;
varepsilon = 10^-5;
n=0;
while max(max(abs(u-un)))>varepsilon
    u=un;
    v=vn;
    u(3:ny-2,1)=10*1/sqrt(2); %L 
    u(3:nx-2,ny)=0;%R
    u(1,3:nx-2)=0; %b
    u(ny,3:nx-2)=0;%t
    v(3:ny-2,1)=0*1/sqrt(2); %L
    v(3:nx-2,ny)=0;%R
    v(1,3:nx-2)=0;%b
    v(ny,3:nx-2)=0;%t
    
    for i=2:(ny-1)
        for j=2:(nx-1)
        un(i,j)=u(i,j)- (dt/2/dx) * u(i,j)*(u(i+1,j) -u(i-1,j)) - ...
            (dt/2/dy) * v(i,j)*(u(i,j+1)-u(i,j-1)) + ...
            (nu*dt/dx^2) *(u(i+1,j)-2*u(i,j)+u(i-1,j)) + ...
            (nu*dt/dy^2) * (u(i,j+1)-2*u(i,j)+u(i,j-1));
        vn(i,j)=v(i,j)- (dt/2/dx) * u(i,j)*(v(i+1,j) -v(i-1,j)) - ...
            (dt/2/dy) * v(i,j)*(v(i,j+1)-v(i,j-1)) + ...
            (nu*dt/dx^2) *( v(i+1,j)-2*v(i,j)+v(i-1,j)) + ...
            (nu*dt/dy^2) * (v(i,j+1)-2*v(i,j)+v(i,j-1));
        end
    end
    un(1:ny,1)=u(1:ny,1);
    un(1,1:nx)=u(1,1:nx);
    un(ny,1:nx)=u(ny,1:nx);
    un(1:nx,ny)=u(1:nx,ny);
    vn(1:ny,1)=v(1:ny,1);
    vn(1,1:nx)=v(1,1:nx);
    vn(ny,1:nx)=v(ny,1:nx);
    vn(1:nx,ny)=v(1:nx,ny);
    if mod(n,500)==0
        tiledlayout(2,4)
        nexttile([2,2])
        image(flip(sqrt(un.^2+vn.^2)),"CDataMapping",'scaled');
        colorbar; caxis([0 1])
        hold off;
        nexttile([2,2])
        quiver(x,y,v,u); 
        hold off;

%     surf(x,y,sqrt(un.^2+vn.^2))
        pause(0.001)
    end
    n=n+1;
end