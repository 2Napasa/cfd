clc;clear;clf; close all; format long;
n = 101;
dx = 1/(n-1);
dt = 10^-4;
alpha = 3;
u0 = 0;
ux = 0;
m = 1;
yps = 10^-5;
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
%     next(n-1)=0;
    for i = n-1:-1:1  
        next(i)=next(i+1)*a(i+1)+b(i+1);
    end
    prev(1)=0; prev(n)=0;
    curr(1)=0; curr(n)=0;
    next(1)=0; next(n)=0;
    length(curr)
    length(prev)
    length(next)
    if mod(ctr,200)==0 
        plot(x,next,'b');
        axis([0 1 -8 8]);
        pause(0.1);
    end
    prev  = curr;
    curr = next;
    m = max(abs(prev-curr));
    next = curr+1;
    ctr = ctr+1;
end

function y = ut0(x)
    y = sin(2*pi*x)+5*sin(3*pi*x);
    %y = 50*sin(pi*x);
end

function y = upt0(x)
    y = 3*sin(5*pi*x);
    %y = 0;
end