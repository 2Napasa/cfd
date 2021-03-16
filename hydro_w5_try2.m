clc;clear;clf; close all; format long;
n = 201;
iter = 0;
dx = 1/(n-1);
dt = 0.5*dx*dx;
eps = 0.00001;
alpha = 9.0;
f=zeros([1 n]);
x=0:dx:1;
oldP = ut0(x);
oldP1=oldP-dt*upt0(x);
newP=oldP+1;
oldP(1)=0; oldP(n)=0;
oldP1(1)=0; oldP1(n)=0;
newP(1)=0; newP(n)=0;
m=1;
ctr = 0;
while m>eps
    oldP(1)=0; oldP(n)=0;
    oldP1(1)=0; oldP1(n)=0;
    newP(1)=0; newP(n)=0;
    for i = 2:n-1
        newP(i)=2*oldP(i)-oldP1(i)+alpha*dt*dt*((oldP(i+1)-2*oldP(i)+oldP(i-1))/(dx*dx)+f(i));
    end
    
    m=max(abs(oldP-oldP1));
    if mod(ctr,1000)==0
        plot(x,newP); axis([0 1 -8 8]);hold off; pause(0.01);
    end
    oldP1=oldP;
    oldP=newP;
    ctr=ctr+1;
end


function y = ut0(x)
    y = sin(2*pi*x)+5*sin(3*pi*x);
end

function y = upt0(x)
    y = 3*sin(5*pi*x);
end