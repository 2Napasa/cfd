% Rauan Kelesbekov
% Given Burgers 2D system of PDE's 
% u_t + u * u_x + v * u_y = (1/Re) * (u_xx + u_yy)
% v_t + u * v_x + v * v_y = (1/Re) * (v_xx + v_yy)
% To be solved using Fractional Step Method
%
% Below written is the function FSM that takes Domain size (Lx,Ly), 
% number of nodes (nx,ny), timestep (dt), Reynold's number (Re), 
% tolerance (varepsilon), inlet (lefthole), outlet(righthole).
% 
% The main goal of this code is to use it in solving Navier-Stokes
% equation, as calling the premade function is more p efficient than
% writing the Burgers Solver from scratch. 
%
% Main issues produced here were the implementation of mixed boundary
% conditions, i.e. Neumann/Dirichlet.
% 
% Moreover, while solving in each direction using Thomas algorithm, the
% zero values on boundaries were recalculated were some values lower than 
% a*10^-6 due to in-built matlab tolerance, therefore the local 
% tolerance TOL = 1e-6 was introduced to make the
% neglectibel? values that are lower than 10 to the -6 equal to zero. 
% During debugging process, new TOL showed itself to work only on
% boundaries, thus not affecting the solution itself. 
% 
% The test on inlet/outlet boundaries could be seen in this video
% https://www.youtube.com/watch?v=VLDY8GTCSmc&feature=youtu.be
% 
% below is the function. 

function [u,v] = FSM(Lx,Ly,nx,ny,dt,Re,varepsilon,lefthole,righthole)
    TOL = 1e-6;
    dx = Lx / (nx-1);
    dy = Ly / (ny-1);
    % initialise FDS
    [x, y] = meshgrid(0:dx:Lx,0:dx:Ly);
    %initialize initial
    us       = zeros(nx,ny);
    vs       = zeros(nx,ny);
    u       = zeros(nx,ny);
    v       = zeros(nx,ny);
    % boundary
    us(:,1)  = 0;
    us(:,ny) = 0;
    us(1,lefthole)  = 1;
    % u(nx,righthole) = u(nx-1,righthole);
    vs(:,1)  = 0;
    vs(:,ny) = 0;
    vs(1,:)  = 0;
    vs(nx,:) = 0;
    %initialize intermediate
    uplus05 = us;
    uplus1  = us;
    vplus05 = vs;
    vplus1  = vs;
    %
    t = 0;
    m = 1;
    ctr  = 0;
    while m > varepsilon
        %initialise current timestep calc U
        usplus05 = us;
        %thomas by x direction
        for j=2:ny-1
            a = zeros(1,nx);
            b = zeros(1,nx);
            a(2)=0; b(2)=us(1,j); %boundary conditions left
            % forward thomas
            for i = 2:nx-1
                Ax = -1/2/Re/dx^2 + us(i,j)/4/dx;
                Bx = 1/dt+1/Re/dx^2;
                Cx = -1/2/Re/dx/dx - us(i,j)/4/dx;
                Dx = (us(i,j)/dt...
                    +(1/2/Re/dx/dx) * (us(i+1,j)-2*us(i,j)+us(i-1,j))...
                    -(vs(i,j)/4/dx)  * (us(i+1,j) + us(i-1,j))...
                    +(1/Re/dy/dy)   * (us(i,j+1)-2*us(i,j)+us(i,j-1))...
                    -(vs(i,j)/2/dy)  * (us(i,j+1) + us(i,j-1)));
                a(i+1) = (-Ax) / (Bx + Cx * a(i));
                b(i+1) = (Dx - Cx * b(i)) / (Bx + Cx *a(i));
            end
            % recalc bdry for current 
            if ismember(j,righthole) %boundary conditions left
                uplus05(nx,j) = b(nx)/(1-a(nx));
            else
                uplus05(nx,j) = 0; %us(nx,j);
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
            a(2)=0; b(2)=us(i,1); % bot boundary all Dirichlet
            % forward thomas
            for j = 2:ny-1
                Ax = -1/2/Re/dy/dy + vs(i,j)/4/dy;
                Bx = 1/dt+1/Re/dy/dy;
                Cx = -1/2/Re/dy/dy - vs(i,j)/4/dy;
                Dx = ((-1/2/Re/dy/dy) * (us(i,j+1)-2*us(i,j)+us(i,j-1))...
                    +(vs(i,j)/4/dy) * (us(i,j+1) + us(i,j-1))...
                    +uplus05(i,j)/dt);
                a(j+1) = (-Ax) / (Bx + Cx * a(j));
                b(j+1) = (Dx - Cx * b(j)) / (Bx + Cx *a(j));
    %             plot(a.','r-'); hold on; plot(b.','b-');  drawnow;
            end

            % recalc bdry for current 
    %         uplus1(nx,righthole) = 0;%b(nx)/(1-a(nx));
            uplus1(:,1)  = 0; % top boundary all Dirichlet
            for j = ny-1:-1:1     
                uplus1(i,j)=uplus1(i,j+1)*a(j+1)+b(j+1);
    %             quiver(x,y,(uplus1).',(v).'); hold on; 
    % %             arrow = zeros(nx,ny);
    %             arrow(i,j) = 1;
    %             quiver(x,y,arrow.',arrow*0); hold off;
    %             axis([-0.1 1.1 -0.1 1.1]);
    %             drawnow; 
            end 
        end
        m = max(max(uplus1-us));
        us = uplus1;
        
            % SOLVING FOR V % % SOLVING FOR V % % SOLVING FOR V % % SOLVING FOR V % % SOLVING FOR V % 
        vsplus05 = v; 
        vsplus1 = v;
        for j = 2:ny-1
            a = zeros(1,nx);
            b = zeros(1,nx);
            % FORWARD THOMAS
            % v(1,:) = v(2,:) * a2 + b2,
            % v(1,:) = array => a2 = 0; b2 = array = v(1,:);
            a(2) = 0; b(2) = v(1,j);
            for i=2:nx-1%CENTRAL DIFFERENCE USED dv/dx = [v(i+1) - v(i-1)]/2/dx
                A = -1/2/Re/dx^2 + u(i,j)/4/dx;
                B = 1/dt+1/Re/dx^2;
                C = -1/2/Re/dx/dx - u(i,j)/4/dx;
                D = (v(i,j)/dt...
                    +(1/2/Re/dx/dx) * (v(i+1,j)-2*v(i,j)+v(i-1,j))...
                    -(v(i,j)/4/dx)  * (v(i+1,j) + v(i-1,j))...
                    +(1/Re/dy/dy)   * (v(i,j+1)-2*v(i,j)+v(i,j-1))...
                    -(v(i,j)/2/dy)  * (v(i,j+1) + v(i,j-1)));
                a(i+1) = (-A) / (B + C * a(i));
                b(i+1) = (D - C * b(i)) / (B + C * a(i));
            end
            % backward THOMAS
            % if WALL => USE DirichLET
            % v(ny) = 0; WALL EVERYWHERE
            vsplus05(nx,j) = v(nx,j);
            for i = nx-1:-1:1
                vsplus05(i,j) = vsplus05(i+1,j) * a(i+1) + b(i+1);
            end
        end
        for i = 2:nx-1
            a = zeros(1,ny);
            b = zeros(1,ny);
            % v(1,:) = v(2,:) * a2 + b2,
            % v(1,:) = array => a2 = 0; b2 = array = v(1,:);
            a(2) = 0; b(2) = v(i,1); % USING LOWER WALL
            for j=2:ny-1
                A = -1/2/Re/dy/dy + v(i,j)/4/dy;
                B = 1/dt+1/Re/dy/dy;
                C = -1/2/Re/dy/dy - v(i,j)/4/dy;
                D = (-1/2/Re/dy/dy) * (v(i,j+1)-2*v(i,j)+v(i,j-1))...
                    +(v(i,j)/4/dy) * (v(i,j+1) + v(i,j-1))...
                    +vsplus05(i,j)/dt;
                a(j+1) = (-A) / (B + C * a(j));
                b(j+1) = (D - C * b(j)) / (B + C *a(j));
            end
            %backward thomas
            % if WALL => USE DirichLET
            % v(ny) = 0; WALL EVERYWHERE
            vsplus1(i,1)  = 0; % USING UPPER WALL
            for j = nx-1:-1:1
                vsplus1(i,j) = vsplus1(i,j+1)*a(j+1)+b(j+1);
            end
        end    
        v = vsplus1;
        
        % plot when difference is enough
%         if mod(ctr,every)==0
    %         quiver(x,y,(us).',(vs).'); hold off;
    %         surf(y,x,flipud(u));
    %         image(u.',"CDataMapping",'scaled'); colorbar;
    %         fprintf('max delta u = %.8f\n',m);
%             drawnow;
%         end    
        % reassigning + increment ctr and time
        t = t + dt; 
        ctr = ctr + 1;
    end
    u = us;
    v = vs;
    u(u<0 & u>-TOL) = 0;
    v(v<0 & v>-TOL) = 0;
    u(u>0 & u<TOL) = 0;
    v(v>0 & v<TOL) = 0;
    
end
%{
% clc; clear; close;
% format long;
% figure('Renderer', 'painters', 'Position', [900 100 1400 700]);
% Lx = 1;
% Ly = 1;
% nx = 81;
% ny = 81;
% dx = Lx / (nx-1);
% dy = Ly / (ny-1);
% dt = 0.01; 
% Re = 400;
% varepsilon = 10^-5;
% % initialise FDS
% x = 0:dx:Lx;
% y = 0:dy:Ly;
% [x, y] = meshgrid(0:dx:Lx,0:dx:Ly);
% Lx = 1; Ly = 1; nx = 81; ny = 81; dt = 0.01; Re = 400; varepsilon = 1e-5; 
% lefthole = floor(3*(ny-1)/4):ny-2;
% righthole = 2:floor((ny-1)/4);
% righthole = floor(3*(ny-1)/4):ny-2;
% [u, v] = FSM(Lx,Ly,nx,ny,dt,Re,varepsilon,lefthole,righthole);
% quiver(x,y,u',v')
% surf(x,y,flipud(u))
%}