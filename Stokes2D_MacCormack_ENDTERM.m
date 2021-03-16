% Rauan 17BD110477
% given 2dimensional NavierStokes system of five equations for compressible
% fluid:
% 1) rho_t + (rhou)_x + (rhov)_y = 0                     - Continuity equation
% 2) (rhou)_t + (rhouu + p -tauxx)_x + (rhouv - tauyx)_y = 0      - x Momentum 
% 3) (rhov)_t + (rhouv - tauxy)_x + (rhovv + p -tauyy)_y = 0      - y Momentum
% 4) E_t + [(E + p)u + qx - utauxx - vtauxy]_x 
%   + [(E+p)v + qy - utauyx - vtauyy]_y = 0                 - Energy equation
% 5) p = rhoRT                                          - Ideal gas law
% 
% auxiliary equations, constants:
% 1.    _ denotes partial derivative in following variable direction
% 2.    E = rho(e + V^2/2) energy conservation
% 3     mu = mu0 (T/T0)^1.5 * (T0+110)/(T+110) Sutherland's Law viscosity 
% 4.    lambda = -2/3 mu    second viscosity due to Stokes hypothesis
% 5.    tauxx = lambda(nablaV) + 2 mu u_x
%       tauyy = lambda(nablaV) + 2 mu v_y
%       tauxy = tauyx = mu (u_y + v_x)
% 6.    rho - density
% 7.    u - velocity in x direction 
% 8.    v - velocity in y direction
% 9.    Pr = mu*cp/k = 0.71         Prandtl number
% 10.   k - thermal conductivity
% 11.   M = 4                       Mach number used for initial velocity
% 12.   qx = -k T_x; qy = -k T_y    heat flux vectors
% 13.   p - pressure
% 14.   R - gas constant 
% 15.   e = cv * T from eq (2.) 
% 16.   cv = R / (gamma - 1)        Specific volume heat constant
% 17.   gamma = 1.4                 Specific heat ratio for perfect air
% 18.   cp = cv * gamma             Specific pressure heat constant
% 19.   V = sqrt(u^2 + v^2)         Absolute velocity 
% 20.   T = 288                     Kelvins           
% 21.   SS = 340                    Speed of sound m/s
% 
% boundary conditoins
%                   
%               T = Tw; u = Mach * soundSpeed; v = 0; P = Patm 
%           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           %                                               
%           %
% T = Tw    %                                                   T_x = 0 
% u = M*SS  %                                                   u_x = 0
% v = 0     %                                                   v_x = 0
% P = Patm  %                                                   P_x = 0
%           %
%           %
%           %
%           %
%           %
%           %
%           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               T = Tw; u = 0; v = 0; P_y = 0
% 
% Method:
% The initial system of N-S equations to be solved using MacCormack 2 step
% predictor-corrector method:
% 
% Consider vectors:
% U = [rho; rhou; rhov; E]
% E = [rhou; rhouu + p - tauxx; rhouv - tauxy; (E+p)u - utauxx - vtauxy + qx]
% F = [rhov; rhouv - tauyx; rhovv + p - tauyy; (E+p)v - utauyx - vtauyy + qy]
%
% The following MacCormack method implies
% 1) Ubar = U - dt/dx (E_{i+1,j} - E_ij) - dt/dy (F_{i,j+1} - F(i,j))
% calculating Ebar and Fbar similarly to E F and rhobar, ubar, vbar, pbar, 
% Tbar, mubar, lambdabar, kbar intermediate values to be used in corrector
% step below
% 2) Unew = 0.5(Ubar + U - dt/dx(Ebar_ij - Ebar_{i-1.j}) - dt/dy(Fbar_ij - Fbar_{i,j-1}))
% extracting corrected rho, u, v, p, T, mu, lambda, k from Unew
% repeating the procedure above
% 
% Results/Conclusion:
% The code below used predetermined time domain of 1000 iterations and
% grid of 71x71 nodes
% Code itself generates a GIF file named 'stokes2d_MacCormack.gif' with 20
% frames
% 
% Results display the Temperature, Density, Velocity and Mach number for
% the given domain. As an experiment, different approach to right outflow
% boundary condition was tested - extrapolation, with three points 
% considered - (un+un-2)=2*un-1, by expressing un alternative boundary
% condition to be used, however the results did not show drastic change,
% therefore it was decided to stay with Neuman conditions. 

% clear cache, close figures, clear command line, 
clc; close all; clear;

% time
ntBreak = 1e3;  % infinite loop break threshhold
ctr     = 1;    % step counter
time    = 0;    % time in seconds

%  gif creator
figure('Renderer', 'painters', 'Position', [1000 200 1600 500]);%plotcoord
filename = 'stokes2d_MacCormack.gif';
FramesNu = 20;
plotRate = ntBreak/FramesNu; % plot only 20 frames

% nodes 
nx = 71; 
ny = 71; 

% stream conditions
M_stream    = 4;           % sea level Mach number
a_stream    = 340;         % sound speed
p_stream    = 101325;      % p atmospheric
T_stream    = 288;         % temp of stream
wallLen     = 1e-5;        
streamHeight= wallLen;     

% constants

mu_0    = 1.78e-5;  % viscosity 
T_0     = 288;      % initial temp
Pr      = 0.71;     % Prandtl nu. for perfect air
R       = 287;      % Specific gas constant 

% gas constants
gamma       = 1.4;              % Specific heat ratio for perfect air
c_v         = R / (gamma - 1);  % Specific volume heat constant 
c_p         = c_v * gamma;      % Specific pressure heat constant

% initial density from p = rhoRT
rho_stream  = p_stream / (R * T_stream);

% step lengths
dx = wallLen/(nx - 1);
dy = streamHeight/(ny - 1);
x = 0:dx:wallLen; 
y = 0:dy:streamHeight; 
delta_t = 5e-11;

% temp, rho and press
p       = ones(ny,nx) * p_stream;
rho     = ones(ny,nx) * rho_stream;
T       = ones(ny,nx) * T_stream;
T(1,:)  = T_stream;  % initial boundaries
T(:,nx) = T_stream;
T(ny,:) = T_stream;

% velocities
u       = ones(ny,nx) * M_stream * a_stream;
u(1,:)  = 0; % bot boundary at wall
v       = zeros(ny,nx);

% viscosity, second viscosity
mu      = viscosity(T, mu_0, T_0);
lambda  = - 2/3 * viscosity(T, mu_0, T_0); % second viscosity https://www.mdpi.com/2076-3417/9/24/5444/htm
k       = TempConductivity(Pr, c_p, mu);

% predictor variables
U1_pred = zeros(ny, nx);
U2_pred = zeros(ny, nx);
U3_pred = zeros(ny, nx);
U5_pred = zeros(ny, nx); % U4 reserved for z coodrinate
rho_pred = zeros(ny, nx);
u_pred = zeros(ny, nx);
v_pred = zeros(ny, nx);
T_pred = zeros(ny, nx);
p_pred = zeros(ny, nx);

% time loop
while ctr <= ntBreak 
    U1 = rho; 
    U2 = rho .* u; 
    U3 = rho .* v; 
    U5 = rho .* (c_v * T + (u.^2 + v.^2)/2);  
    
    % Predictor
    [E1, E2, E3, E5] = E(rho, u, p, v, T, mu, lambda, k, c_v, dx, dy, 'EPred');
    [F1, F2, F3, F5] = F(rho, u, p, v, T, mu, lambda, k, c_v, dx, dy, 'FPred');
    for i = 2:nx-1
        for j = 2:ny-1
            U1_pred(j,i) = U1(j,i) ...
                - delta_t * ((E1(j,i+1) - E1(j,i))/dx + (F1(j+1,i) - F1(j,i))/dy);
            U2_pred(j,i) = U2(j,i) ...
                - delta_t * ((E2(j,i+1) - E2(j,i))/dx + (F2(j+1,i) - F2(j,i))/dy);
            U3_pred(j,i) = U3(j,i) ...
                - delta_t * ((E3(j,i+1) - E3(j,i))/dx + (F3(j+1,i) - F3(j,i))/dy);
            U5_pred(j,i) = U5(j,i) ...
                - delta_t * ((E5(j,i+1) - E5(j,i))/dx + (F5(j+1,i) - F5(j,i))/dy);
        end
    end
    
    % pred interior from UPred
    [rho_pred(2:ny-1,2:nx-1), u_pred(2:ny-1,2:nx-1), v_pred(2:ny-1,2:nx-1), T_pred(2:ny-1,2:nx-1)] ... 
                = rhouvT(U1_pred(2:ny-1,2:nx-1), ...
                    U2_pred(2:ny-1,2:nx-1), ...
                    U3_pred(2:ny-1,2:nx-1), ...
                    U5_pred(2:ny-1,2:nx-1), ...
                    c_v);
    p_pred(2:ny-1,2:nx-1) = rho_pred(2:ny-1,2:nx-1) * R .* T_pred(2:ny-1,2:nx-1);
    
    % pred Boundary
    [rho_pred, u_pred, v_pred, p_pred, T_pred] = ...
            bndr(rho_pred, u_pred, v_pred, p_pred, T_pred, rho_stream, ...
                M_stream * a_stream, p_stream, T_stream, R);
    
    % pred coeffs
    mu_pred     = viscosity(T_pred, mu_0, T_0);
    lambda_pred = - 2/3 * mu_pred;
    k_pred      = TempConductivity(Pr, c_p, mu_pred);
    
    
    % Corrector
    [E1_p, E2_p, E3_p, E5_p] = ...
                    E(rho_pred, u_pred, p_pred, v_pred, T_pred, mu_pred, ...
                    lambda_pred, k_pred, c_v, dx, dy, 'ECorr');
    [F1_p, F2_p, F3_p, F5_p] = ...
                    F(rho_pred, u_pred, p_pred, v_pred, T_pred, mu_pred, ...
                    lambda_pred, k_pred, c_v, dx, dy, 'FCorr');
    for i = 2:nx-1
        for j = 2:ny-1
            U1(j,i) = 1/2*(U1(j,i) + U1_pred(j,i) - delta_t*((E1_p(j,i) ...
                - E1_p(j,i-1))/dx + (F1_p(j,i) - F1_p(j-1,i))/dy));
            U2(j,i) = 1/2*(U2(j,i) + U2_pred(j,i) - delta_t*((E2_p(j,i) ...
                - E2_p(j,i-1))/dx + (F2_p(j,i) - F2_p(j-1,i))/dy));
            U3(j,i) = 1/2*(U3(j,i) + U3_pred(j,i) - delta_t*((E3_p(j,i) ...
                - E3_p(j,i-1))/dx + (F3_p(j,i) - F3_p(j-1,i))/dy));
            U5(j,i) = 1/2*(U5(j,i) + U5_pred(j,i) - delta_t*((E5_p(j,i) ...
                - E5_p(j,i-1))/dx + (F5_p(j,i) - F5_p(j-1,i))/dy));
        end
    end
    
    % corr interior from Unew
    [rho(2:ny-1,2:nx-1), u(2:ny-1,2:nx-1), v(2:ny-1,2:nx-1), T(2:ny-1,2:nx-1)] ...
                = rhouvT(U1(2:ny-1,2:nx-1), ...
                    U2(2:ny-1,2:nx-1), ...
                    U3(2:ny-1,2:nx-1), ...
                    U5(2:ny-1,2:nx-1), ...
                    c_v);
    p(2:ny-1,2:nx-1) = rho(2:ny-1,2:nx-1)*R.*T(2:ny-1,2:nx-1);
    
    % corr Boundary
    [rho,u,v,p,T] = bndr(rho,u,v,p,T,rho_stream,M_stream*a_stream,p_stream,T_stream,R);
    
    % corr coeffs
    mu = viscosity(T, mu_0, T_0);
    lambda = - 2/3*mu;
    k = TempConductivity(Pr, c_p, mu);
     
    % creating GIF with number of FramesNu images
    if mod(ctr,plotRate) == 1
        % plotter
        tiledlayout(1,3);
        nexttile; 
        contourf(x,y,T,10,'r');  colorbar; hold on;
        axis([-0.1e-5 1.1e-5 -0.1e-5 1.1e-5])
        title('Temperature [K]');
        nexttile; 
        contourf(x,y,rho,10,'b');  colorbar; hold off;
        axis([-0.1e-5 1.1e-5 -0.1e-5 1.1e-5])
        title('Density [kg/m3]');
        nexttile; 
        M = sqrt(u.^2 + v.^2)./sqrt(gamma*R*T);
        contourf(x,y,M,10,'b'); colorbar; hold on;
        axis([-0.1e-5 1.1e-5 -0.1e-5 1.1e-5])
%         title('Mach Number');
%         nexttile; 
        norm = u.^2 + v.^2;
        quiver(x,y,u,v,.2,'k'); colorbar; hold off; 
        axis([-0.1e-5 1.1e-5 -0.1e-5 1.1e-5])
        title('velocity field [Mach]');
        drawnow;
            % gif creator
            frame = getframe(1);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if ctr == 1
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
              imwrite(imind,cm,filename,'gif','WriteMode','append');
            end
    end
    ctr = ctr + 1;
    time = time + delta_t;
end

% some useful functions:

function [E1, E2, E3, E5] = E(rho, u, p, v, T, mu, lambda, k, c_v, dx, dy, PredCorrFlag)
    E1 = rho .* u;
    
    tau_xx = Tauxx(u, v, lambda, mu, dx, dy, PredCorrFlag);
    E2 = rho.*u.^2 + p - tau_xx;
    
    tau_xy = Tauxy(u, v, mu, dx, dy, PredCorrFlag);
    E3 = rho.*u.*v - tau_xy;
    
    qx = qxFun(T, k, dx, PredCorrFlag);
    E5 = (rho.*(c_v*T + (u.^2 + v.^2)/2) + p).*u - u.*tau_xx - v.*tau_xy + qx;
end

function [F1, F2, F3, F5] = F(rho, u, p, v, T, mu, lambda, k, c_v, dx, dy, PredCorrFlag)
    F1 = rho.*v;
    
    tau_yx = Tauxy(u, v, mu, dx, dy, PredCorrFlag);
    F2 = rho.*u.*v - tau_yx;
    
    tau_yy = Tauyy(u, v, lambda, mu, dx, dy, PredCorrFlag);
    F3 = rho.*v.^2 + p - tau_yy;
    
    qy = qyFun(T, k, dy, PredCorrFlag);
    F5 = (rho.*(c_v*T + (u.^2 + v.^2)/2) + p).*v - u.*tau_yx - v.*tau_yy + qy;
end

function [rho, u, v, T] = rhouvT(U1, U2, U3, U5, c_v)
    rho = U1;
    u = U2 ./ U1;
    v = U3 ./ U1;
    T = (U5 ./ U1 - ((U2 ./ U1).^2 + (U3 ./ U1).^2) / 2) / c_v;
end

function k = TempConductivity(Pr, c_p, mu)
    k = mu*c_p/Pr;
end

% in functions of stresses below, derivatives at boundaries have only one possible
% difference from forward or rearward to not exceed array indicies
function tau_xx = Tauxx(u, v, lambda, mu, dx, dy, PredCorrFlag)
    [ny, nx] = size(u);
    u_x = zeros(ny, nx);
    v_y = zeros(ny, nx);
    
    % du/dx:
    if (strcmp(PredCorrFlag, 'EPred'))
        for i = 2:nx
            for j = 1:ny
                u_x(j,i) = (u(j,i) - u(j,i-1)) / dx; 
            end
        end
        u_x(:,1) = (u(:,2) - u(:,1))/dx;
    elseif (strcmp(PredCorrFlag, 'ECorr'))
        for i = 1:nx-1
            for j = 1:ny
                u_x(j,i) = (u(j,i+1) - u(j,i))/dx;
            end
        end
        u_x(:,nx) = (u(:,nx) - u(:,nx-1))/dx;
    else
        error('wrong PredCorrFlag, please enter one of Epred Ecorr.\n')
    end
    
    % dv/dy:
    for i = 1:nx
        for j = 2:ny-1
            v_y(j,i) = (v(j+1,i) - v(j-1,i))/(2*dy);
        end
    end
    v_y(1,:) = (v(2,:) - v(1,:))/dy; 
    v_y(ny,:) = (v(ny,:) - v(ny-1,:))/dy;

    tau_xx = lambda .* (u_x + v_y) + 2 * mu .* u_x;
end

function tau_yy = Tauyy(u, v, lambda, mu, dx, dy, PredCorrFlag)
    [ny, nx] = size(v);
    u_x = zeros(ny, nx);
    v_y = zeros(ny, nx);
    
    % dv/dy:
    if (strcmp(PredCorrFlag, 'FPred'))
        for i = 1:nx
            for j = 2:ny
                v_y(j,i) = (v(j,i) - v(j-1,i))/dy;
            end
        end
        v_y(1,:) = (v(2,:) - v(1,:))/dy; 
    elseif (strcmp(PredCorrFlag, 'FCorr'))
        for i = 1:nx
            for j = 1:ny-1
                v_y(j,i) = (v(j+1,i) - v(j,i))/dy;
            end
        end
        v_y(ny,:) = (v(ny,:) - v(ny-1,:))/dy; 
    else
        error('wrong PredCorrFlag, please enter one of Fpred Fcorr.\n')
    end

    % du/dx
    for i = 2:nx-1
        for j = 1:ny
            u_x(j,i) = (u(j,i+1) - u(j,i-1))/(2*dx);
        end
    end
    u_x(:,1) = (u(:,2) - u(:,1))/dx;
    u_x(:,nx) = (u(:,nx) - u(:,nx-1))/dx;
    
    tau_yy = lambda .* (u_x + v_y) + 2 * mu .* v_y;
end

function tau_xy = Tauxy(u, v, mu, dx, dy, PredCorrFlag)

    [ny, nx] = size(u);
    u_y = zeros(ny, nx);
    v_x = zeros(ny, nx);
    if (strcmp(PredCorrFlag, 'EPred') || strcmp(PredCorrFlag, 'ECorr'))        
        for i = 1:nx
            for j = 2:ny-1
                u_y(j,i) = (u(j+1,i) - u(j-1,i))/(2*dy);
            end
        end
        u_y(1,:) = (u(2,:) - u(1,:))/dy; 
        u_y(ny,:) = (u(ny,:) - u(ny-1,:))/dy; 
        if (strcmp(PredCorrFlag, 'EPred'))
            for i = 2:nx
                for j = 1:ny
                    v_x(j,i) = (v(j,i) - v(j,i-1))/dx;
                end
            end
            v_x(:,1) = (v(:,2) - v(:,1))/dx;
        else
            for i = 1:nx-1
                for j = 1:ny
                    v_x(j,i) = (v(j,i+1) - v(j,i))/dx;
                end
            end
            v_x(:,nx) = (v(:,nx) - v(:,nx-1))/dx;
        end
    elseif (strcmp(PredCorrFlag, 'FPred') || strcmp(PredCorrFlag, 'FCorr'))
        for i = 2:nx-1
            for j = 1:ny
                v_x(j,i) = (v(j,i+1) - v(j,i-1))/(2*dx);
            end
        end
        v_x(:,1) = (v(:,2) - v(:,1))/dx;
        v_x(:,nx) = (v(:,nx) - v(:,nx-1))/dx;
        if (strcmp(PredCorrFlag, 'FPred'))
            for i = 1:nx
                for j = 2:ny
                    u_y(j,i) = (u(j,i) - u(j-1,i))/dy;
                end
            end
            u_y(1,:) = (u(2,:) - u(1,:))/dy;
        else
            for i = 1:nx
                for j = 1:ny-1
                    u_y(j,i) = (u(j+1,i) - u(j,i))/dy;
                end
            end
            u_y(ny,:) = (u(ny,:) - u(ny-1,:))/dy;
        end
    else
        error('wrong PredCorrFlag, please enter one of Epred Ecorr Fpred Fcorr.\n')
    end

    tau_xy = mu .* (u_y + v_x);
end



function qy = qyFun(T, k, dy, PredCorrFlag)
    [ny, nx] = size(T);
    T_y = zeros(ny, nx);
    if (strcmp(PredCorrFlag, 'FPred'))
        for i = 1:nx
            for j = 2:ny
                T_y(j,i) = (T(j,i) - T(j-1,i))/dy;
            end
        end
        T_y(1,:) = (T(2,:) - T(1,:))/dy; 
    elseif (strcmp(PredCorrFlag, 'FCorr'))
        for i = 1:nx
            for j = 1:ny-1
                T_y(j,i) = (T(j+1,i) - T(j,i))/dy; 
            end
        end
        T_y(ny,:) = (T(ny,:) - T(ny-1,:))/dy; 
    else
        error('wrong PredCorrFlag, please enter one of Fpred Fcorr.\n')
    end

    qy = -k .* T_y;
end

function qx = qxFun(T, k, dx, PredCorrFlag)
    [ny, nx] = size(T);
    T_x = zeros(ny, nx);
    if (strcmp(PredCorrFlag, 'EPred'))
        for i = 2:nx
            for j = 1:ny
                T_x(j,i) = (T(j,i) - T(j,i-1))/dx;
            end
        end
        T_x(:,1) = (T(:,2) - T(:,1))/dx;
    elseif (strcmp(PredCorrFlag, 'ECorr'))
        for i = 1:nx-1
            for j = 1:ny
                T_x(j,i) = (T(j,i+1) - T(j,i))/dx;
            end
        end
        T_x(:,nx) = (T(:,nx) - T(:,nx-1))/dx;
    else
        error('wrong PredCorrFlag, please enter one of Epred Ecorr.\n')
    end

    qx = -k .* T_x;
end

function mu = viscosity(T, mu_0, T_0) % Sutherland's Law
    mu = mu_0 * (T/T_0).^(3/2) * (T_0 + 110) ./ (T + 110);
end


function [rho, u, v, p, T] = bndr(rho, u, v, p, T, rho_stream, ...
                                    u_stream, p_stream, T_stream, R)
    % consider boundaries  below, starting at left wall and the following in
    % clockwise direction
    % 
    %
    %           %_%%%%%%%%%%%
    %           %           _
    %           %           %
    %           %           %
    %           %           %
    %           _           _
    %           %%%%%%%%%%%%%
    %
    [ny, nx]        = size(rho);
    % Left
    u(2:ny,1)       = u_stream;
    p(2:ny,1)       = p_stream;
    T(2:ny,1)       = T_stream;
    rho(2:ny,1)     = rho_stream;
    % Upper
    u(ny,2:nx)      = u_stream;
    p(ny,2:nx)      = p_stream;
    T(ny,2:nx)      = T_stream;
    rho(ny,2:nx)    = rho_stream;
    % Right 
    u(2:ny-1,nx)    = u(2:ny-1,nx-1); 
    v(2:ny-1,nx)    = v(2:ny-1,nx-1);
    p(2:ny-1,nx)    = p(2:ny-1,nx-1);
    T(2:ny-1,nx)    = T(2:ny-1,nx-1);
    rho(2:ny-1,nx)  = p(2:ny-1,nx) ./ (R * T(2:ny-1,nx));
    % Bot
    T(1,1:nx)       = T_stream;
    p(1,1:nx)       = p(2,1:nx);
    rho(1,1:nx)     = p(1,1:nx) ./ (R * T(1,1:nx));
end