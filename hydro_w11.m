% Relaxation method. 
% Problem statement
% given Poisson equation u_xx + u_yy = -f_ij
% boundary ctds
%          0
%     %%%%%%%%%%%
%   1 %   u0=0  % 0
%     %         % 
%     %%%%%%%%%%%
%         0
% f_ij = 0
% To be solved using relaxation method
% The method showed itself being way faster 
% than Gauss Seidel method that is having ~6000 iterations 
% vs Relaxation having ~200 iterations on a given domain of 101 vs 101
% nods
% However optimal w has to be found in advance.
% Below is the animation of plots for an optimal w search, it can be seen
% that the converged solution of an u in equation does not change with
% respect to w
% https://www.youtube.com/watch?v=O0hZhU8CSn8

clc;clear;clf; close all;
format long;
every = 1000;
alpha = 1;
radiatortemp = 1;
outtemp = 0;
Lx = 1;
Ly = 1;
nx = 101;
ny = 101;
dx = Lx / (nx-1);
dy = Ly / (ny-1);
dt = (dx^2 + dy^2) * 0.01 ; 
varepsilon = 10^-6;
% initialise general
x = 0:dx:Lx;
y = 0:dy:Ly;
%initialize initial
% theta       = ones(nx,ny)*0;
% 
% % boundary
% theta(1,:) = radiatortemp;
% 
% %initialize intermediate
% thetaplus1  = theta;
% %
% t = 0;
% m = 1;
ctrarr = [];
warr = [];
for w = 1.01:0.01:1.99
    ctr  = 0;
    theta = ones(nx,ny)*0;
    m=1;
    theta(1,:) = radiatortemp;
    thetaplus1 = theta;
while m > varepsilon
%     %initialise current timestep
    theta = thetaplus1;
    thetaplus1(1,:) = 1;
    for i = 2:nx-1
        for j = 2:ny-1
            thetaplus1(i,j) = (w/4) * (theta(i+1,j) + thetaplus1(i-1,j) + theta(i,j+1) + thetaplus1(i,j-1)...
                - 4*(1-1/w)*theta(i,j));
        end
    end
    thetaplus1(1,:) = 1;
%     if mod(ctr,every)==0
%         tic
%         image((theta)','CDataMapping','scaled'); 
%         colormap(jet);
%         caxis([0 1]);colorbar; 
% %         surf(thetaplus1);
%         pause(0.01);
% %         fprintf('iter = %d, delta = %f, %f\n',ctr,m,toc)
% 
%     end
    m = max(max(abs(thetaplus1(3:nx-2,3:ny-2)) - abs(theta(3:nx-2,3:ny-2))));
    ctr = ctr + 1;
end
ctrarr = [ctrarr ctr];
warr = [warr w];
tiledlayout(2,1)
nexttile
plot(warr,ctrarr,'x'); axis([1 2 0 6000])
nexttile
image((thetaplus1)','CDataMapping','scaled');  colorbar; colormap(jet); 
%         surf(thetaplus1);
pause(0.001);
fprintf('iter = %d; w = %f\n',ctr,w)
end
function c = boundary(nx,ny)
    c = 1;
end