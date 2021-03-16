%%%variable declarations
clear; close all; clc;
format long;
figure('Renderer', 'painters', 'Position', [900 100 1400 700]);
nx = 111;
ny = 111;
Lx = 1; Ly = 1;
Re = 100;
dx = Lx/(nx-1); dy = Ly/(ny-1);
dt = 0.1;
[x, y] = meshgrid(0:dx:Lx,0:dx:Ly);
u = zeros(ny,nx); 
v = zeros(ny,nx);
varepsilon = 10^-5;
ctr = 0;

lefthole = 2:ny-1; %floor(3*(ny-1)/4):ny-2;
righthole = 2:ny-1; %2:floor((ny-1)/4);
u(1,lefthole)=1;                    % inlet
u(nx,righthole)=u(nx-1,righthole);  % outlet
unew = u;
vnew = v;
u = zeros(1,nx);
u(1)=1;
for k = 1:100
%     u(1,lefthole)=1;                    % inlet
%     u(nx,righthole)=u(nx-1,righthole);  % outlet
    uplus1 = u;
    a = zeros(1,nx); b = a; 
    a(2) = 0;
    b(2) = 1;
    for i=2:ny-1
%         A = -1/dx/dx/Re;
%         B = 1/dt + 2/dx/dx/Re + u(i)/dx;
%         C = -1/dx/dx/Re - u(i)/dx;
%         D = u(i-1)/dt;
        A = u(i)/2/dx - 1/dx/dx/Re;
        B = 1/dt + 2/Re/dx/dx;
        C = -1/dx/dx/Re - u(i)/2/dx;
        D = u(i)/dt;
        a(i+1) = (-A) / (B + C * a(i));
        b(i+1) = (D - C * b(i)) / (B + C * a(i));
    end
%     plot(a,'r'); hold on; plot(b,'b');hold off; drawnow; pause(0.5)
    uplus1(nx) = b(nx)/(1-a(nx));
%     b(nx)/(1-a(nx))
    for i = nx-1:-1:1
       uplus1(i) = uplus1(i+1)*a(i+1) + b(i+1); 
    end
    plot(0:dx:1,uplus1,'k'); axis([0 1 0 1]); drawnow;
    u = uplus1;
    
    if mod(ctr,100)==0
%         quiver(x,y,u.',v.',0.5,'r'); drawnow;
    end
    ctr = ctr + 1;
end