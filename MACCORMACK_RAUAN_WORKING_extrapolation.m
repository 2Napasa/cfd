% clear cache, close figures, clear command line, 
clc; close all; clear;

% time
ntBreak = 1e3;  % infinite loop break threshhold
ctr     = 1;    % step counter
time    = 0;    % time in seconds

%  gif creator
figure('Renderer', 'painters', 'Position', [1000 200 1200 1200]);%plotcoord
filename = 'stokes2d_MacCormack.gif';
FramesNu = 20;
plotRate = ntBreak/FramesNu; % plot only 20 frames

% nodes 
nx = 71; 
ny = 71; 

% stream conditions
M_stream    = 4;           % sea level Mach number
a_stream    = 340.28;         % sound speed
p_stream    = 101325;      % p atmospheric
T_stream    = 288.16;         % temp of stream
wallLen     = 1e-5;        
streamHeight= wallLen;     

% constants
gamma   = 1.4;      % Specific heat ration for perfect air
mu_0    = 1.78e-5;  % viscosity 
T_0     = 288;      % initial temp
Pr      = 0.71;     % Prandtl nu. for perfect air
R       = 287;      % Specific gas constant 

% gas constants
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
u(1,:)  = 0; % bot boundary
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
        tiledlayout(2,2);
        nexttile; 
        contourf(x,y,T,10,'k');  colorbar;
        title('Temperature');
        nexttile; 
        contourf(x,y,rho,10,'k');  colorbar;
        title('density');
        nexttile; 
        M = sqrt(u.^2 + v.^2)./sqrt(gamma*R*T);
        contourf(x,y,M,10,'b'); colorbar;
        title('Mach Number');
        nexttile; 
        norm = u.^2 + v.^2;
        quiver(x,y,u,v,.4,'r'); colorbar;
        axis([-0.1e-5 1.1e-5 -0.1e-5 1.1e-5])
        title('velocity field');
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


function [rho, u, v, p, T] = bndr(rho, u, v, p, T, rho_stream, u_stream, p_stream, T_stream, R)
    % consider boundaries  below, starting at (1,1) and the following in
    % clockwise direction
    % right side considered with three point extrapolation (un+un-2)=2*un-1
    %
    %           %_%%%%%%%%%%%
    %           %           _
    %           %           %
    %           %           %
    %           %           %
    %           _           _
    %           %_%%%%%%%%%%%
    %
    
    [ny, nx]        = size(rho);
    % Corner
    T(1,1)          = T_stream;
    p(1,1)          = p_stream;
    rho(1,1)        = rho_stream;
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
    u(2:ny-1,nx)    = 2 * u(2:ny-1,nx-1) - u(2:ny-1,nx-2); 
    v(2:ny-1,nx)    = 2 * v(2:ny-1,nx-1) - v(2:ny-1,nx-2);
    p(2:ny-1,nx)    = 2 * p(2:ny-1,nx-1) - p(2:ny-1,nx-2);
    T(2:ny-1,nx)    = 2 * T(2:ny-1,nx-1) - T(2:ny-1,nx-2);
    rho(2:ny-1,nx)  = p(2:ny-1,nx) ./ (R * T(2:ny-1,nx));
    % Bot
    T(1,2:nx)       = T_stream;
    p(1,2:nx)       = 2 * p(2,2:nx) - p(3,2:nx);
    rho(1,2:nx)     = p(1,2:nx) ./ (R * T(1,2:nx));
end