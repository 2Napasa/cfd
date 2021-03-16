% Rauan Kelesbekov
% Chem proc - reaction of two reagents producing heat
% Given the problem
% A + B = C + D + T, leads to system of PDE's
% T_t = k*A*B + alpha^2 * T_xx - u_T * T_x
% A_t + u_A * A_x = -kAB
% B_t + u_B * B_x = -kAB
% and the parallel reaction 
% D_t + u_A * D_x = -k1DE
% E_t + u_B * D_x = -k1DE
% F_t             =  k1DE
% k1 = k1 * gamma^(T); is changing due to Arrhenius equation
%  
% A,B constantly supplied on left (1/4) and right (3/4) of domain
% vertically, Likewise D and E. 
% after discretizing came up to numerical solution below.
% k = 300 to speed up the reaction; alpha = 0.01 to lower the diffusion. 
% u_A = 3 to the right side -> step backward used for CFL stability
% u_B = -3 to the left side -> step forward  used for CFL stability
% of the Implicit method. 
%
% The final result displays the reagents A and B coming together by
% advection and heat T diffusion propagation after being produced after 
% A and B combined. 
% 
% Temperature has advecion part to transport faster to the neighbouring 
% process to speed up the reaction
%
% Had to play with coefficients to speed up / slow down reaction D + E = F
% to make visualizing clearer
% 
% the result could be seen in the file attached
%
% the code: 

clc; clear; close;
nodes  =  101;
L  =  1;
k  =  300;
alpha  =  0.01;
dx = L / (nodes - 1);
dt = 10^ - 6;
Ac = zeros(1,nodes);
Bc = zeros(1,nodes);
Cc = zeros(1,nodes);
Dc = zeros(1,nodes);
Ec = zeros(1,nodes);
Fc = zeros(1,nodes);
Tc = zeros(1,nodes);

An = Ac;
Bn = Bc;
Cn = Cc;
Dn = Dc;
En = Ec;
Fn = Fc;
Tn = Tc; 
varepsilon  =  10^ - 6;
every  =  1600;
ctr  =  0;
x = 0:dx:L;
leftA = floor(1 * nodes/8);
rightB = floor(3 * nodes/8);
leftD = floor(5 * nodes/8);
rightE = floor(7 * nodes/8);
Bc(rightB) = 1; % IC for B reagent
Ac(leftA) = 1;
m  =  1;
ua = 3;%3; 
ub = -3;%-3;
h = figure;
axis tight manual
filename = 'testChem.gif';
k1 = ones(1,nodes)*k/2;
gamma = 1.1;
% c3h8 + 5*o2 = 3*co2 + 4*h2o + t
for n = 1:650000
    %constant supply
    for i = 2:nodes-1 %change in reagents
        An(i) = dt*(-k*1*Ac(i)^1*Bc(i)^1+(-ua)*(Ac(i)-Ac(i-1))/(dx)) + Ac(i);
        Bn(i) = dt*(-k*1*Ac(i)^1*Bc(i)^1+(-ub)*(Bc(i+1)-Bc(i))/(dx)) + Bc(i);
%         Cn(i) = dt*(k*3)*Ac(i)^1*Bc(i)^5 + Cc(i);
%         Dn(i) = dt*(k*4)*Ac(i)^1*Bc(i)^5 + Dc(i);
        Tn(i) = dt*((k*160)*(Ac(i)^1*Bc(i))^5+(alpha/dx^2)*(Tc(i+1)-2*Tc(i)+Tc(i-1))+ (-ua/2)*(Tc(i)-Tc(i-1))/(dx)) + Tc(i);
        Dn(i) = dt*(-k1(i)*1*Dc(i)^1*Ec(i)^1+(-ua)*(Dc(i)-Dc(i-1))/(dx)) + Dc(i);
        En(i) = dt*(-k1(i)*1*Dc(i)^1*Ec(i)^1+(-ub)*(Ec(i+1)-Ec(i))/(dx)) + Ec(i);
        Fn(i) = dt*(k1(i)*1/40*Dc(i)^1*Ec(i)^1) + Fc(i);
        k1(i) = k1(i) * gamma^((Tn(i))/400);
    end
    
    An(rightB)=An(rightB-1);
    Bn(leftA)=Bn(leftA+1);
    Dn(rightE)=Dn(rightE-1);
    En(leftD)=En(leftD+1);
%     Cn(nodes)=Cn(nodes-1);
%     Dn(nodes)=Dn(nodes-1);
    Tn(nodes)=0;
    if mod(ctr,every) == 0 %plot faster
        plot(x,An,'r',x,Bn,'r'); hold on; 
        plot(x,Dn,'b',x,En,'b'); hold on;
        plot(x,Fn,'k'); hold on;
        plot(x,Tn,'g'); hold off;
        hold off; 
        axis([dx L - dx 0 1]); 
        drawnow;   
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if n == 1
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
        else
          imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
        
    end 

    Ac = An;
    Ac(leftA) = 1;
    Bc = Bn;
    Bc(rightB) = 1; % adding water

    Dc = Dn;
    Dc(leftD) = 1;
    Ec = En;
    Ec(rightE) = 1;

    Fc = Fn;
    Tc = Tn;
    ctr = ctr + 1;
    
end