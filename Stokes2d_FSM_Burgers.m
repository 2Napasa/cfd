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
nx = 25;
ny = 25;
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
    us=u;
    vs=v;
    usplus05 = u;
    usplus1  = u;
    vsplus05 = v;
    vsplus1  = v;
    %{
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
    %}
    m = 1;
ctrBurgers  = 0;
    while m > varepsilon
        %initialise current timestep calc U
        us = usplus1;
        %thomas by dx direction
        for j=2:ny-1
            a = zeros(1,nx);
            b = zeros(1,nx);
            a(2)=0; b(2)=us(1,j); %boundary conditions bottom
            % forward thomas
            for i = 2:nx-1
                uxplus05    = (us(i+1,j) + us(i,j))/2;
                uxminus05   = (us(i,j) + us(i-1,j))/2;
                uyplus05    = (us(i,j+1) + us(i,j))/2;
                uyminus05   = (us(i,j) + us(i,j-1))/2;
                Ax = -1/2/Re/dx^2 + us(i,j)/4/dx;
                Bx = 1/dt+1/Re/dx^2;
                Cx = -1/2/Re/dx/dx-us(i,j)/4/dx;
                Dx = (us(i,j)/dt...
                    +(1/2/Re/dx/dx) * (us(i+1,j)-2*us(i,j)+us(i-1,j))...
                    -(vs(i,j)/2/dx) * (uxplus05 - uxminus05)...
                    +(1/Re/dy/dy) * (us(i,j+1)-2*us(i,j)+us(i,j-1))...
                    -(vs(i,j)/1/dy) * (uyplus05 - uyminus05));
                a(i+1) = (-Ax) / (Bx + Cx * a(i));
                b(i+1) = (Dx - Cx * b(i)) / (Bx + Cx *a(i));
            end
            % recalc bdry for current 
            usplus1(:,1)  = 0;
            usplus1(:,ny) = 0;

            % backward thomas
            for i = nx-1:-1:1      
                usplus1(i,j)=usplus1(i+1,j)*a(i+1)+b(i+1);
    %             quiver(x,y,(uplus1./sqrt(uplus1.^2 + vplus1.^2).*sqrt(uplus1.^2 + vplus1.^2))',(v./sqrt(uplus1.^2 + vplus1.^2).*sqrt(uplus1.^2 + vplus1.^2))'); drawnow;
            end 
        end

        %solving in j direction
        usplus05 = usplus1;
        %thomas by dy direction
        for i=2:nx-1
            a = zeros(1,ny);
            b = zeros(1,ny);
            a(2)=0; b(2)=us(i,1); % left boundary all Direchlet
            % forward thomas
            for j = 2:ny-1
                usplus05    = (us(i,j+1) + us(i,j))/2;
                usminus05   = (us(i,j) + us(i,j-1))/2;
    % % % % % % % % % % %             vsplus05    = (vs(i,j+1) + vs(i,j))/2;
    % % % % % % % % % % %             vsminus05   = (vs(i,j) + vs(i,j-1))/2;
                Ax = -1/2/Re/dy/dy + vs(i,j)/4/dy;
                Bx = 1/dt+1/Re/dy/dy;
                Cx = -1/2/Re/dy/dy-vs(i,j)/4/dy;
                Dx = ((-1/2/Re/dy/dy) * (us(i,j+1)-2*us(i,j)+us(i,j-1))...
                    +(vs(i,j)/2/dy) * (usplus05-usminus05)...
                    +usplus05(i,j)/dt);
                a(j+1) = (-Ax) / (Bx + Cx * a(j));
                b(j+1) = (Dx - Cx * b(j)) / (Bx + Cx *a(j));
            end
            % recalc bdry for current 
            usplus1(1,lefthole)  = 1;
            usplus1(nx,righthole) = b(i)/(1-a(i));
            % backward thomas
            for j = ny-1:-1:1      
                usplus1(i,j)=usplus1(i,j+1)*a(j+1)+b(j+1);
            end 
        end

        % plot when difference is enough
        if mod(ctrBurgers,every)==0
            quiver(x,y,(us).',(vs).');
            fprintf('max delta u = %.8f\n',m)
            drawnow;
        end
        m = max(max(usplus1-us));



        %initialise current timestep calc V
        vs = vsplus1;
        %thomas by dx direction
        for j=2:ny-1
            a = zeros(1,nx);
            b = zeros(1,nx);
            a(2)=0; b(2)=vs(1,j); %boundary conditions bottom
            % forward thomas
            for i = 2:nx-1
                vxplus05    = (vs(i+1,j) + vs(i,j))/2;
                vxminus05   = (vs(i,j) + vs(i-1,j))/2;
                vyplus05    = (vs(i,j+1) + vs(i,j))/2;
                vyminus05   = (vs(i,j) + vs(i,j-1))/2;
                Ax = -1/2/Re/dx^2 + us(i,j)/4/dx;
                Bx = 1/dt+1/Re/dx^2;
                Cx = -1/2/Re/dx/dx-us(i,j)/4/dx;
                Dx = (vs(i,j)/dt...
                    +(1/2/Re/dx/dx) * (vs(i+1,j)-2*vs(i,j)+vs(i-1,j))...
                    -(vs(i,j)/2/dx) * (vxplus05 - vxminus05)...
                    +(1/Re/dy/dy) * (vs(i,j+1)-2*vs(i,j)+vs(i,j-1))...
                    -(vs(i,j)/1/dy) * (vyplus05 - vyminus05));
                a(i+1) = (-Ax) / (Bx + Cx * a(i));
                b(i+1) = (Dx - Cx * b(i)) / (Bx + Cx *a(i));
            end
            % recalc bdry for current 
            vsplus1(:,1)  = 0;
            vsplus1(:,ny) = 0;

            % backward thomas
            for i = nx-1:-1:1      
                vsplus1(i,j)=vsplus1(i+1,j)*a(i+1)+b(i+1);
    %             quiver(x,y,(uplus1./sqrt(uplus1.^2 + vplus1.^2).*sqrt(uplus1.^2 + vplus1.^2))',(v./sqrt(uplus1.^2 + vplus1.^2).*sqrt(uplus1.^2 + vplus1.^2))'); drawnow;
            end 
        end

        %solving in j direction
        vsplus05 = vsplus1;
        %thomas by dy direction
        for i=2:nx-1
            a = zeros(1,ny);
            b = zeros(1,ny);
            a(2)=0; b(2)=us(i,1); % left boundary all Direchlet
            % forward thomas
            for j = 2:ny-1
                vsplus05    = (vs(i,j+1) + vs(i,j))/2;
                vsminus05   = (vs(i,j) + vs(i,j-1))/2;
    % % % % % % % % % % %             vsplus05    = (vs(i,j+1) + vs(i,j))/2;
    % % % % % % % % % % %             vsminus05   = (vs(i,j) + vs(i,j-1))/2;
                Ax = -1/2/Re/dy/dy + vs(i,j)/4/dy;
                Bx = 1/dt+1/Re/dy/dy;
                Cx = -1/2/Re/dy/dy-vs(i,j)/4/dy;
                Dx = ((-1/2/Re/dy/dy) * (vs(i,j+1)-2*vs(i,j)+vs(i,j-1))...
                    +(vs(i,j)/2/dy) * (vsplus05-vsminus05)...
                    +vsplus05(i,j)/dt);
                a(j+1) = (-Ax) / (Bx + Cx * a(j));
                b(j+1) = (Dx - Cx * b(j)) / (Bx + Cx *a(j));
            end
            % recalc bdry for current 
            vsplus1(1,lefthole)  = 0;
            vsplus1(nx,righthole) = 0;
            % backward thomas
            for j = ny-1:-1:1      
                vsplus1(i,j)=vsplus1(i,j+1)*a(j+1)+b(j+1);
            end 
        end

        % plot when difference is enough
        if mod(ctrBurgers,every)==0
            quiver(x,y,(us).',(vs).');
            fprintf('max delta u = %.8f\n',m)
            drawnow;
        end
        m = max(max(usplus1-us));




        % reassigning + increment ctr and time
        ctrBurgers = ctrBurgers + 1;
    end
    vs = vn; 
    us = un;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       solving Poisson eq                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    m = 1; % pressure convergence parameter
    %boundary conditions
    p(1,:)=p(2,:);  
    p(nx,:)=p(nx-1,:);
    p(:,ny)=p(:,ny-1);
    p(:,1)=p(:,2);
    p(1,lefthole)=1;     % inlet pressure
    p(nx,righthole)=0;         % outlet pressure
    pplus1  = p;
    % calculating free term 
    f = zeros(nx,ny);
    for j = 2:ny-1
        for i = 2:nx-1
            usplus05    = (us(i+1,j) + us(i,j))/2;
            usminus05   = (us(i,j) + us(i-1,j))/2;
            vsplus05    = (vs(i,j+1) + vs(i,j))/2;
            vsminus05   = (vs(i,j) + vs(i,j-1))/2;
            f(i,j) = - (rho/dt) ...
                * ((usplus05-usminus05)/dx ...
                    + (vsplus05-vsminus05)/dy);
        end
    end
    % relaxation method, solving Poisson equation, finding pressure field
    w = 1.95;   % optimal relaxation parameter
    while m > varepsilon
        p = pplus1;
        for i = 2:nx-1
            for j = 2:ny-1
                pplus1(i,j) = (w/4)*(p(i+1,j) + p(i,j+1) ...
                    - 4*(1-1/w)*p(i,j) + dx^2 * f(i,j)...
                    + pplus1(i-1,j) + pplus1(i,j-1));
            end
        end
        m = max(max(abs(pplus1) - abs(p)));
    end
    p = pplus1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    correcting u,v considering pressure                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    uplus1 = us;
    vplus1 = vs;
    for i = 2:nx-1 
        for j = 2:ny-1 
            % staggered mesh
            pxplus005   = (p(i+1,j) + p(i,j))/2;
            pxminus05   = (p(i,j) + p(i-1,j))/2;                        
            pyplus005   = (p(i,j+1) + p(i,j))/2;
            pyminus05   = (p(i,j) + p(i,j-1))/2;
            uplus1(i,j) = - dt / rho * (pxplus005-pxminus05)/dx + us(i,j);
            vplus1(i,j) = - dt / rho * (pyplus005-pyminus05)/dy + vs(i,j);
        end
    end
    delta = max(max(abs(uprev-u)))+max(max(abs(vprev-v)));
    u = uplus1;
    v = vplus1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Fancy plots                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if mod(n,every)==0 
        q = (u.^2 + v.^2); q = sqrt(q);     % vector field length
        unormed = (u./q); vnormed = (v./q); % normalizing vector field
        contourf(x,y,p.',10,'w-'); hold on; % pressure field
        colorbar;                 
        caxis([min(min(p)) max(max(p))]);   % color axis            
        contour(x,y,q.',10,'k-');  hold on; % velocity countour
        quiver(x,y,unormed.',vnormed.',0.4,'k');    % velocity field
        axis([-0.1 1.1 -0.1 1.1]); hold on;
        fprintf('%d. pmax = %f umax = %f vmax = %f delta = %.6f\n',...
            n,max(max(abs(p))),max(max(abs(u))),max(max(abs(v))),delta);
        colorbar vert; 
        
        drawnow;
    end
    n = n + 1;
end
