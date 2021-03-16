% Heat transfer 2D
% Rauan Kelesbekov
% Problem statement 
% Given PDE of heat transfer in thin rod with initial temperature and boundary conditions.
% u_t=alpha*(u_xx + u_yy)
% u(x=0)=1
% u(x=L)=0
% u(y=0)=0
% u(y=L)=0
% u(t=0)=0
% for this problem alpha =1, L = 1 used. 
% FSM method was used. and split into two sets of one dimensional problems
% (by x and by y directions separately)
% The Thomas method is used to solve the above 
% A(i)x(i-1)+B(i)x(i)+C(i)x(i+1)=D(i) based on assumption that 
% x(i)=a(i+1)x(i+1)+b(i+1) in forward sweep where 
% a(i+1)=-B(i)/(A(i)a(i)+C(i))
% b(i+1)=(F(i)-A(i)b(i))/(A(i)a(i)+C(i)) and i=2,3,...,n  found in backward sweep
% fixing the y coordinates and performing the above
% The initial assumption is that a(2)=0 b(2)=1 from X boundary conditions
% using theta(i) = theta(i+1)a(i+1) + b(i+1)
% i.e. left/right boundary conditions. and using FSM for coefficients ABCD
% calculating intermediate funciton theta
% then fixing the x coordinate and performing similar process with
% different initial conditions (found from previous) and  boundary for
% upper and lower
% 
% The method showed itself to be absoultely stable, however intermediate
% theta(n+1/2) is required to calculate the u(n+1). 
% 
% below is the code:
clc; clear; close;
format long;
figure('Renderer', 'painters', 'Position', [900 100 1400 700]);
every = 100;
Lx = 1;
Ly = 1;
nx = 11;
ny = 11;
dx = Lx / (nx-1);
dy = Ly / (ny-1);
dt = 0.01; 
alpha = 1;
Re = 100;
varepsilon = 10^-5;
% initialise FDS
x = 0:dx:Lx;
y = 0:dy:Ly;
[x, y] = meshgrid(0:dx:Lx,0:dx:Ly);
%initialize initial
u       = zeros(nx,ny);
v       = zeros(nx,ny);
% boundary
lefthole = floor(3*(ny-1)/4):ny-2;
righthole = 2:floor((ny-1)/4);
lefthole = 2:nx-1; righthole = 2:nx-1;
% righthole = 2:ny;
u(:,1)  = 0;
u(:,ny) = 0;
u(1,lefthole)  = 1;
% u(nx,righthole) = u(nx-1,righthole);
v(:,1)  = 0;
v(:,ny) = 0;
v(1,:)  = 0;
v(nx,:) = 0;
%initialize intermediate
uplus05 = u;
uplus1  = u;
vplus05 = v;
vplus1  = v;
%
t = 0;
m = 1;
ctr  = 0;
while m > varepsilon
    %initialise current timestep calc U
    usplus05 = u;
    %thomas by x direction
    for j=2:ny-1
        a = zeros(1,nx);
        b = zeros(1,nx);
        a(2)=0; b(2)=u(1,j); %boundary conditions left
        % forward thomas
        for i = 2:nx-1
            Ax = -1/2/Re/dx^2 + u(i,j)/4/dx;
            Bx = 1/dt+1/Re/dx^2;
            Cx = -1/2/Re/dx/dx - u(i,j)/4/dx;
            Dx = (u(i,j)/dt...
                +(1/2/Re/dx/dx) * (u(i+1,j)-2*u(i,j)+u(i-1,j))...
                -(v(i,j)/4/dx)  * (u(i+1,j) + u(i-1,j))...
                +(1/Re/dy/dy)   * (u(i,j+1)-2*u(i,j)+u(i,j-1))...
                -(v(i,j)/2/dy)  * (u(i,j+1) + u(i,j-1)));
            a(i+1) = (-Ax) / (Bx + Cx * a(i));
            b(i+1) = (Dx - Cx * b(i)) / (Bx + Cx *a(i));
        end
        % recalc bdry for current 
        if ismember(j,righthole) %boundary conditions left
            uplus05(nx,j) = b(nx)/(1-a(nx));
        else
            uplus05(nx,j) = u(nx,j);
        end
        % backward thomas
        for i = nx-1:-1:1      
            uplus05(i,j)=uplus05(i+1,j)*a(i+1)+b(i+1);
        end 
        
    end
    
    %solving in y direction
    %thomas by dy direction
    for i=2:nx
        a = zeros(1,ny);
        b = zeros(1,ny);
        a(2)=0; b(2)=u(i,1); % bot boundary all Direchlet
        % forward thomas
        for j = 2:ny-1
            Ax = -1/2/Re/dy/dy + v(i,j)/4/dy;
            Bx = 1/dt+1/Re/dy/dy;
            Cx = -1/2/Re/dy/dy - v(i,j)/4/dy;
            Dx = ((-1/2/Re/dy/dy) * (u(i,j+1)-2*u(i,j)+u(i,j-1))...
                +(v(i,j)/4/dy) * (u(i,j+1) + u(i,j-1))...
                +uplus05(i,j)/dt);
            a(j+1) = (-Ax) / (Bx + Cx * a(j));
            b(j+1) = (Dx - Cx * b(j)) / (Bx + Cx *a(j));
%             plot(a.','r-'); hold on; plot(b.','b-');  drawnow;
        end
        
        % recalc bdry for current 
        uplus1(:,1)  = 0; % top boundary all Direchlet
        for j = ny-1:-1:1     
            uplus1(i,j)=uplus1(i,j+1)*a(j+1)+b(j+1);
        end 
    end
    m = max(max(uplus1-u));
    u = uplus1;
    % plot when difference is enough
    if mod(ctr,every)==0
        quiver(x,y,(u).',(v).'); hold off;
        surf(-x,y,u.');
        fprintf('max delta u = %.8f\n',m);
        drawnow;
    end
    
   
    
    
    % reassigning + increment ctr and time
    t = t + dt; 
    ctr = ctr + 1;
end
% quiver(x,y,(u).',(v).');
% fprintf('max delta u = %.8f\n',m)