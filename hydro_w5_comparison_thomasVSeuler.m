clc;clear;clf; close all; format long;
n = 101;
dx = 1/(n-1);
dt = 10^-4;
alpha = 3;
u0 = 0;
ux = 0;
m = 1;
yps = 10^-8;
x = 0:dx:1;
prev = ut0(x);
curr = upt0(x)*dt+prev;
next = prev + 1;
a = curr;
b = curr;
ctr = 0;
A = (-1*alpha^2)/(dx^2);
B = (2*alpha^2/dx^2)+(1/dt^2);
C = -1*alpha^2/(dx^2);
plotctr=0;
set(gcf,'WindowState','fullscreen')
t=0;
while m > yps
    prev(1)=0; prev(n)=0;
    curr(1)=0; curr(n)=0;
    next(1)=0; next(n)=0;
    a(2) = 0; b(2) = u0;
    D = (2*curr-prev)/(dt^2);
    %forw
    for i = 2:n-1
        a(i+1) = (-A)/(B + C*a(i));
        b(i+1) = (D(i)-C*b(i))/(B+C*a(i));
    end
    %backw
    next(n)=0;
    next(n-1)=0;
    for i = n-1:-1:1  
        next(i)=next(i+1)*a(i+1)+b(i+1);
    end
    prev(1)=0; prev(n)=0;
    curr(1)=0; curr(n)=0;
    next(1)=0; next(n)=0;
    newP=prev;
    for i = 2:n-1
        newP(i)=2*curr(i)-prev(i)+alpha*dt*dt*((curr(i+1)-2*curr(i)+curr(i-1))/(dx*dx));
    end
    if plotctr>100
        ana = sin(2*pi*x)*cos(6*pi*t)+5*sin(3*pi*x)*cos(9*pi*t)+sin(t*pi*x)*sin(15*pi*t)/5/pi;
        tiledlayout(3,2)
        nexttile
        plot(x,next,'b'); hold off;
        title('Implicit')
        axis([0 1 -8 8]);
        nexttile()
        plot(x,next-ana); hold off;
        title('Difference implicit')
        axis([0 1 -1 1]);
        nexttile
        plot(x,newP,'k'); hold off;
        title('Explicit')
        axis([0 1 -8 8]);
        nexttile()
        plot(x,newP-ana); hold off;
        title('Difference explicit')
        axis([0 1 -1 1]);
        nexttile
        plot(x,ana,'r'); hold off;
        title('Analytic')
        axis([0 1 -8 8]);
        pause(0.01);
        plotctr=0;
    end
    prev  = curr;
    curr = next;
    m = max(abs(prev-curr));
    next = curr+1;
    ctr = ctr+1;
    plotctr = plotctr+m;
    t=t+dt;
end

function y = ut0(x)
    y = sin(2*pi*x)+5*sin(3*pi*x);
    %y = 50*sin(pi*x);
end

function y = upt0(x)
    y = 3*sin(5*pi*x);
    %y = 0;
end