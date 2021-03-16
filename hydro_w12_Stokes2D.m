% u_t + u * u_x + v * u_y = -1/rho * p_x + u_xx + u_yy
% v_t + u * v_x + v * v_y = -1/rho * p_y + v_xx + v_yy
% u_x + v_y = 0

clc;clear;clf; close all;
format long;
figure('Renderer', 'painters', 'Position', [10 10 900 900]);
every = 20;
inletspeed = 1;
Lx = 1; % length
Ly = 1; % length
nx = 41; % nu of nodes by x
ny = 41; % nu of nodes by y
dx = Lx / (nx-1); 
dy = Ly / (ny-1);
dt = dx^4; 
varepsilon = 10^-4;
Re = 2;
rho =  0.1;
% initialise general
% x = 0:dx:Lx;
% y = 0:dy:Ly;
[x,y] = meshgrid(0:dx:Lx,0:dx:Ly);
u = zeros(nx,ny);
v = zeros(nx,ny);
p = zeros(nx,ny);

u(nx,2:(ny-1)*0.25)=0;
u(1,(ny-1)*0.75:ny-1)=inletspeed;
ctr = 1;
us = u;
vs = v; 
for k = 1:100
    % u_t + u * u_x + v * u_y =  + u_xx + u_yy to be solved
%     u(:,ny)=0;
%     u(:,1)=0;
%     v(:,1)=0;
%     v(:,ny)=0;
%     u(nx,:)=u(nx-1,:);
%     u(1,:)=1;
%     v(1,:)=0;
%     v(nx,:)=0;
    
    us = u; vs = v;
    for i = 2:nx-1
        for j = 2:ny-1
%             uxplus05    = (u(i+1,j) + u(i,j))/2;
%             uxminus05   = (u(i,j) + u(i-1,j))/2;
%             vxplus05    = (v(i+1,j) + v(i,j))/2;
%             vxminus05   = (v(i,j) + v(i-1,j))/2;
%             uyplus05    = (u(i,j+1) + u(i,j))/2;
%             uyminus05   = (u(i,j) + u(i,j-1))/2;
%             vyplus05    = (v(i,j+1) + v(i,j))/2;
%             vyminus05   = (v(i,j) + v(i,j-1))/2;
%             us(i,j) = u(i,j) ...
%                 - (dt/dx) * u(i,j) * (uxplus05 - uxminus05) ...
%                 - (dt/dy) * v(i,j) * (uyplus05 - uyminus05) ...
%                 + (1/Re*dt/dx^2) * (u(i+1,j) - 2*u(i,j) + u(i-1,j))...
%                 + (1/Re*dt/dy^2) * (u(i,j+1) - 2*u(i,j) + u(i,j-1));
%             vs(i,j) = v(i,j)...
%                 - (dt/dx) * u(i,j) * (vxplus05 - vxminus05) ...
%                 - (dt/dy) * v(i,j) * (vyplus05 - vyminus05) ...
%                 + (1/Re * dt/dx^2) * (v(i+1,j) - 2*v(i,j) + v(i-1,j))...
%                 + (1/Re * dt/dy^2) * (v(i,j+1) - 2*v(i,j) + v(i,j-1));

%                   ????? ? ?????????? ?????? ??????
                    us(i,j)= u(i,j) ...
                        - dt ...
                        * (max(u(i,j),0) * (u(i,j) - u(i-1,j))/dx ...
                        + min(u(i,j),0) * (u(i+1,j) - u(i,j))/dx)...
                        - dt ...
                        * (max(v(i,j),0) * (u(i,j) - u(i,j-1))/dy ...
                        + min(v(i,j),0) * (u(i,j+1) - u(i,j))/dx) ... 
                        + (1/Re*dt/dx^2) * (u(i+1,j)-2*u(i,j)+u(i-1,j)) ...
                        + (1/Re*dt/dy^2) * (u(i,j+1)-2*u(i,j)+u(i,j-1));
                    vs(i,j)= v(i,j) ... 
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

%     nexttile;
    
%     image((us),'CDataMapping','scaled'); colorbar; 
    
    m = 1;
    p(1,:)=p(2,:);
    p(nx,:)=p(nx-1,:);
    p(:,ny)=p(:,ny-1);
    p(:,1)=p(:,2);
    p(1,(ny-1)*0.75:ny-1)=1; %inlet
    p(nx,2:(ny-1)*0.25)=0;  %outlet
    pplus1  = p;
    w = 1.95;       %relaxation parameter
    f = zeros(nx,ny);
    for j = 2:ny-1
        for i = 2:nx-1
            usplus05    = (us(i+1,j) + us(i,j))/2;
            usminus05   = (us(i,j) + us(i-1,j))/2;
            vsplus05    = (vs(i,j+1) + vs(i,j))/2;
            vsminus05   = (vs(i,j) + vs(i,j-1))/2;
            f(i,j) = - (rho/dt) * ((usplus05-usminus05)/dx + (vsplus05-vsminus05)/dy); %????
        end
    end
    while m > varepsilon
    %     %initialise current timestep
        p = pplus1;
        for i = 2:nx-1
            for j = 2:ny-1
                pplus1(i,j) = (w/4)*(p(i+1,j) + p(i,j+1) - 4*(1-1/w)*p(i,j) + dx^2 * f(i,j)...
                    + pplus1(i-1,j) + pplus1(i,j-1));
            end
        end
        m = max(max(abs(pplus1) - abs(p)));
    end
    p = pplus1;
    
    uplus1 = u;
    vplus1 = v;
    for i = 2:nx-1
        for j = 2:ny-1
            pxplus05    = (p(i+1,j) + p(i,j))/2;
            pxminus05   = (p(i,j) + p(i-1,j))/2;                        
            pyplus05    = (p(i,j+1) + p(i,j))/2;
            pyminus05   = (p(i,j) + p(i,j-1))/2;
            uplus1(i,j) = -dt / rho * (pxplus05-pxminus05)/dx + us(i,j);
            vplus1(i,j) = -dt / rho * (pyplus05-pyminus05)/dy + vs(i,j);
            
        end
        
    end
    
    u = uplus1;
    v = vplus1;
    u(nx,2:(ny-1)*0.25)=u(nx-1,2:(ny-1)*0.25); %outlet
    u(1,(ny-1)*0.75:ny-1)=inletspeed;    %inlet
 
    
    
    if mod(ctr,every)==0
        contourf(x,y,p.',10,'w-'); colorbar; axis([-0.1 1.1 -0.1 1.1]); hold on;
        drawnow;
%         quiver(x,y,u.',v.');    
        quiver(x,y,(u.'),(v.'));  axis([-0.1 1.1 -0.1 1.1]);
        drawnow;
        fprintf('%d. max p = %f max u = %f\n',ctr,max(max(p)),max(max(abs(u))));
        
    end
    if max(max(p)) > 300090000
            break;
    end

    ctr = ctr + 1;
end