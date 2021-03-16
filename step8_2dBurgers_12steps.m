%%%variable declarations
clear all
nx = 41;
ny = 41;
nt = 5000;
every = nt / 20;
c=1;
dx = 2.0/(nx-1);
dy = 2.0/(ny-1);
sigma = .009;
Re = 1; 
dt = dx^4;

x = linspace(0,2,nx);
y = linspace(0,2,ny);

u = zeros(ny,nx); %%create a 1xn vector of 1's
v = zeros(ny,nx);
un = zeros(ny,nx);
vn = zeros(ny,nx);


%%%Assign initial conditions

% % u(.5/dy:1/dy+1,.5/dx:1/dx+1)=2; %%set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2
% % v(.5/dy:1/dy+1,.5/dx:1/dx+1)=2; %%set hat function I.C. : u(.5<=x<=1 && .5<=y<=1 ) is 2

%%%Plot Initial Condition
 %%the figsize parameter can be used to produce different sized images
                  
[X, Y] = meshgrid(x,y);                            


inletspeed = 1;
ctr = 1; 
for n=1:nt+1
    un=u;
    vn=v;
    for i=2:(ny-1)
        for j=2:(nx-1)
        u(i,j)=u(i,j) ...
            - (dt/dx) * u(i,j) * (u(i,j) -u(i-1,j)) ...
            - (dt/dy) * v(i,j) * (u(i,j)-u(i,j-1)) ...
            + (1/Re*dt/dx^2) * (u(i+1,j)-2*u(i,j)+u(i-1,j)) ...
            + (1/Re*dt/dy^2) * (u(i,j+1)-2*u(i,j)+u(i,j-1));
        end
    end
    for i = 2:nx-1
        for j = 2:ny-1
            v(i,j)=v(i,j) ...
                - (dt/dx) * u(i,j) * (v(i,j) -v(i-1,j)) ...
                - (dt/dy) * v(i,j) * (v(i,j)-v(i,j-1)) ...
                + (1/Re*dt/dx^2) * (v(i+1,j)-2*v(i,j)+v(i-1,j)) ...
                + (1/Re*dt/dy^2) * (v(i,j+1)-2*v(i,j)+v(i,j-1));
       end
   end
   u(nx,2:(ny-1)*0.25)=inletspeed;
   u(1,(ny-1)*0.75:ny-1)=inletspeed;
   if mod(ctr,every)==0
       quiver(x,y,u.',v.'); 
       drawnow;
   end
   ctr = ctr + 1;
end

