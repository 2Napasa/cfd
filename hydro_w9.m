clc;clear;clf; close all;
format long;
every = 100;
alpha = 1;
radiatortemp = 90;
outtemp = 0;
Lx = 1;
Ly = 1;
nx = 101;
ny = 101;
dx = Lx / (nx-1);
dy = Ly / (ny-1);
dt = (dx^2 + dy^2) * 0.1 ; 
varepsilon = 10^-6;
% initialise general
x = 0:dx:Lx;
y = 0:dy:Ly;
%initialize initial
theta       = ones(nx,ny)*20;

% boundary
theta(:,1:2)  = outtemp;                          %out  top
theta(:,2) = (theta(:,1) + theta(:,3))/2;   %wall top
theta(:,ny-1:ny) = outtemp;                            %out  bot
theta(:,ny-1) = (theta(:,ny) + theta(:,ny-2))/2; %wall bot
theta(1:2,:) = outtemp;                           %out  left
theta(2,:) = (theta(1,:) + theta(3,:))/2;   %wall left
theta(1:2,floor(1*nx/5):floor(2*nx/5))  = radiatortemp; %left radiator1
theta(1:2,floor(3*nx/5):floor(4*nx/5))  = radiatortemp; %left radiator2
theta(nx-1:nx,:)=outtemp;                                 %out  right
theta(nx-1,:) = (theta(nx,:) + theta(nx-2,:))/2;    %wall right
theta(nx,floor(1*nx/5):floor(2*nx/5)) = ...
         theta(nx-2,floor(nx/5):floor(2*nx/5)); %right hole1
theta(nx,floor(3*nx/5):floor(4*nx/5)) = ...
         theta(nx-2,floor(3*nx/5):floor(4*nx/5));   %right hole2
theta(nx-1,floor(1*nx/5):floor(2*nx/5)) = ...
         theta(nx-2,floor(nx/5):floor(2*nx/5)); %right hole1
theta(nx-1,floor(3*nx/5):floor(4*nx/5)) = ...
         theta(nx-2,floor(3*nx/5):floor(4*nx/5));   %right hole2

%initialize intermediate
thetaplus05 = theta;
thetaplus1  = theta;
%
t = 0;
m = 1;
ctr  = 0;
while m > varepsilon
%     %initialise current timestep
    theta = thetaplus1;
    thetaplus1(:,1:2)  = outtemp;                          %out  top
    thetaplus1(:,2) = (thetaplus1(:,1) + thetaplus1(:,3))/2;   %wall top
    thetaplus1(:,ny-1:ny) = outtemp;                            %out  bot
    thetaplus1(:,ny-1) = (thetaplus1(:,ny) + thetaplus1(:,ny-2))/2; %wall bot
    thetaplus1(1:2,:) = outtemp;                           %out  left
    thetaplus1(2,:) = (thetaplus1(1,:) + thetaplus1(3,:))/2;   %wall left
    thetaplus1(1:2,floor(1*nx/5):floor(2*nx/5))  = radiatortemp; %left radiator1
    thetaplus1(1:2,floor(3*nx/5):floor(4*nx/5))  = radiatortemp; %left radiator2
    thetaplus1(nx:nx,:)=outtemp;                                 %out  right
    thetaplus1(nx-1,:) = (thetaplus1(nx,:) + thetaplus1(nx-2,:))/2;    %wall right
    thetaplus1(nx-1,floor(1*nx/5):floor(2*nx/5)) = ...
             thetaplus1(nx-3,floor(nx/5):floor(2*nx/5)); %right hole1
    thetaplus1(nx-1,floor(3*nx/5):floor(4*nx/5)) = ...
             thetaplus1(nx-3,floor(3*nx/5):floor(4*nx/5));   %right hole2
    thetaplus1(nx,floor(1*nx/5):floor(2*nx/5)) = ...
             thetaplus1(nx-3,floor(nx/5):floor(2*nx/5)); %right hole1
    thetaplus1(nx,floor(3*nx/5):floor(4*nx/5)) = ...
             thetaplus1(nx-3,floor(3*nx/5):floor(4*nx/5));   %right hole2
    for i = 2:nx-1
        for j = 2:ny-1
            thetaplus1(i,j) = ...
                + (alpha ^ 2 * dt / dx^2) * (theta(i+1,j) - 2*theta(i,j) + theta(i-1,j))...
                + (alpha ^ 2 * dt / dy^2) * (theta(i,j+1) - 2*theta(i,j) + theta(i,j-1)) ...
                + theta(i,j);  ...
%                 - theta(i,j)*0.1* (theta(i+1,j) - theta(i-1,j))/2/dx*dt; %% convection term
        end
    end

    if mod(ctr,every)==0
        tic
        image((thetaplus1)','CDataMapping','scaled'); 
        colormap(jet);
        caxis([outtemp radiatortemp]);colorbar; 
%         surf(thetaplus1);
        pause(1);
        fprintf('iter = %d, delta = %f, %f\n',ctr,m,toc)

    end
    m = max(max(abs(thetaplus1(3:nx-2,3:ny-2)) - abs(theta(3:nx-2,3:ny-2))));
    ctr = ctr + 1;
end
image((thetaplus1)','CDataMapping','scaled');  colorbar; 
%         surf(thetaplus1);
pause(0.1);
fprintf('iter = %d, delta = %f\n',ctr,m)
function c = boundary(nx,ny)
    c = 1;
end