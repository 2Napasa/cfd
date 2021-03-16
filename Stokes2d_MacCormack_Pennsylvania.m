clear; close all; clc;
format short;
figure('Renderer', 'painters', 'Position', [1000 200 1000 500]);%plotcoord
nx = 21;
ny = 21;
nt = 10;                       % used during testing, no longer req
every = floor(nt / 1 * 10000);  % plot frequency
Lx = 0.00001;
Ly = 5*Lx;
dx = Lx/(nx-1);
dy = Ly/(ny-1);
sigma = .009;
mu=.01;
Re = 1/mu;
dt = dx^3;            % CFL condition
x = linspace(0,1,nx);
y = linspace(0,1,ny);
[x, y] = meshgrid(x,y);  
p = zeros(nx,ny);                   
varepsilon = 1e-7;
delta = 1;          % to enter the main loop
figure(1);
plotctr = 0;

% initial and boundary cond's
T = ones(nx,ny)*288;
R = 8.31;
T0=288;
% pressures
patm = 1e5;
p(1,:) = patm;      % pressure left dirichlet
p(:,ny) = patm;     % pressure top dirichlet
p(nx,:) = 0;        % pressure out dirichlet
p(:,1) = p(:,2);    % pressure bot normal

Cv = (3/2) * R; 
Cp = (5/2) * R;
e = T ./ Cv;
rho0 = patm/R/T0;
rho = ones(nx,ny) * rho0;                           % not needed?
u = zeros(ny,nx);                              % gives error cant divide by zero
v = zeros(ny,nx);
uleft = 1; utop = 1; ubot = 0;
vleft = 0; vtop = 0; vbot = 0;
u(1,:)  = uleft;    v(1,:)  = vleft;
u(:,ny) = utop;     v(:,ny) = vtop;
u(:,1)  = ubot;     v(:,1)  = vbot;




rhonew = rho; 
rhounew = rho;
rhovnew = rho;
pnew = p;
rhou = rho.*u;
rhov = rho.*v;
rhostar = rho;
rhoustar = rhou;
rhovstar = rhov;

t = 0;
ctr = 0; 
c1 = dt/dx; 
c2 = dt/dy;
c3 = mu*dt/dx/dx;
c4 = mu*dt/dy/dy;
c5 = mu*dt/12/dx/dy;
c = R*T;
for time = 1:100
    for i = 2:nx-1
        rhostar = rho; rhoustar = rhou; rhovstar = rhov; 
        for j = 2:ny-1
            rhostar(i,j) = rho(i,j) ...
                -c1 * (rhou(i+1,j) - rhou(i,j))...
                -c2 * (rhov(i,j+1) - rhov(i,j));
            rhoustar(i,j) = rhou(i,j)...
                -c1 * (rhou(i+1,j) * u(i+1,j) + c(i+1,j)^2 * rho(i+1,j) ...
                        - rhou(i,j) * u(i,j) + c(i,j)^2 * rho(i,j))...
                -c2 * (rhou(i,j+1) * v(i,j+1) - rhou(i,j) * v(i,j))...
                +(4*c3/3) * (u(i+1,j) - 2*u(i,j) + u(i-1,j))...
                +c4 * (u(i,j+1) - 2* u(i,j) + u(i,j-1))...
                +c5 * (v(i+1,j+1) + v(i-1,j-1) - v(i+1,j-1)-v(i-1,j+1));
            rhovstar(i,j) = rhov(i,j)...
                -c1 * (rhou(i+1,j) * v(i+1,j) - rhou(i,j) * v(i,j))...
                -c2 * (rhov(i,j+1) * v(i,j+1) + c(i,j+1)^2 * rho(i,j+1)...
                        - rhov(i,j) * v(i,j) + c(i,j) * rho(i,j))...
                +c3 * (v(i+1,j) - 2* v(i,j) + v(i-1,j))...
                +(4*c4/3) * (v(i,j+1) - 2* v(i,j) + v(i,j-1))...
                +c5 * (u(i+1,j+1) - u(i-1,j-1) - u(i+1,j-1) - u(i-1,j+1));
        end
    end
    ustar = rhoustar./rhostar; vstar = rhovstar./rhostar;
    rhoustarSQR = rhoustar .* ustar; rhovstarSQR = rhovstar.*vstar;
    cstar = R.*T;
    rhonew = rho; rhounew = rhou; rhovnew = rhov; 
    
    rhostar(1,:) = rho(1,:) - dt/2/dx * (-rho(3,:) + 4*rho(2,:) .* u(2,:) - 3 * rho(1,:) .* u(1,:));
    rhonew(1,:) = 0.5 * (rho(1,:) + rhostar(1,:) ...
        - dt/2/dx * (rhostar(3,:).*ustar(3,:) ...
        + 4 * rhostar(2,:) .* ustar(2,:) - 3 * rhostar(1,:).* ustar(1,:)));
    for i = 2:nx-1
        for j = 2:ny-1
            rhonew(i,j) = 0.5 * (...
            rho(i,j) + rhostar(i,j)...
            -c1 * (rhoustar(i,j) - rhoustar(i-1,j))...
            -c2 * (rhovstar(i,j) - rhovstar(i,j-1))...
                );
            rhounew(i,j) = 0.5 * (...
                rhou(i,j) + rhoustar(i,j)...
                -c1 * (rhoustarSQR(i,j) + cstar(i,j)^2 * rhostar(i,j)...
                        - rhoustarSQR(i-1,j) + cstar(i-1,j)^2 * rhostar(i-1,j))...
                -c2 * (rhoustar(i,j) * vstar(i,j) - rhoustar(i,j-1) * vstar(i,j-1))...
                +(4*c3/3) * (ustar(i+1,j) - 2*ustar(i,j) + ustar(i-1,j))...
                +c4 * (ustar(i,j+1) - 2*ustar(i,j) + ustar(i,j))...
                +c5 * (vstar(i+1,j+1) + vstar(i-1,j-1) - vstar(i+1,j-1) - vstar(i-1,j+1))...
                );
            rhovnew(i,j) = 0.5 * (...
                rhov(i,j) + rhovstar(i,j)...
                -c1 * (rhoustar(i,j) * vstar(i,j) - rhoustar(i-1,j) * vstar(i,j))...
                -c2 * (rhovstarSQR(i,j) + cstar(i,j)^2 * rhostar(i,j)...
                        - rhovstarSQR(i,j-1) + cstar(i,j-1)^2 * rhostar(i,j-1))...
                +c3 * (vstar(i+1,j) - 2 * vstar(i-1,j))...
                +(4*c4/3) * (vstar(i,j+1) - 2*vstar(i,j) + vstar(i,j-1))...
                +c5 * (ustar(i+1,j+1) + ustar(i-1,j-1) -ustar(i+1,j-1) - ustar(i-1,j+1))...
                );
        end
    end
    
    
    rho = rhonew;
    u = rhounew./rho;
    v = rhovnew./rho;
    u(1,:) = 1; u(nx,:) = u(nx-1,:); u(:,1) = 0; u(:,ny) = 1;
    v(1,:) = 0; v(nx,:) = v(nx-1,:); v(:,1) = 0; v(:,ny) = 0;
    rhou = rho.*u;
    rhov = rho.*v;
    
    
    pnew = rho.*R.*T;
    quiver(rhou.',rhov.'); %hold on;
%     quiver(rhoustar.',rhovstar.');
    drawnow;
    
end