clc;clear all;clf;
L = 1;
n = 100;
dx = L / n;
ld = 3; 
dt = sqrt(dx^2/ld);
x = 0:dx:L;
u_prev = sin(2*pi*x)+5*sin(3*pi*x);
u = dt*(3*sin(5*pi*x)) + u_prev;
u_new = u + 1;
vareps = 10^-8;
delta = 1; 
ctr = 0;
mplot = 1;
while delta > vareps
 	u(1)=0;
    u(n)=0;
    u_prev(1)=0;
    u_prev(n)=0;
	for i=2:n-1
		u_new(i) = ld*dt*dt*(u(i-1)-2*u(i)+u(i+1))/dx/dx + 2*u(i) - u_prev(i);
    end
    u(1)=0;
    u(n)=0;
    u_prev(1)=0;
    u_prev(n)=0;
    if mplot > 9
        plot(x,u); 
        axis([0 1 -50 50])
        pause(0.01);
        mplot=0;
    end
	delta = max(abs(u-u_new));
    u_prev = u;
	u = u_new;
    
    mplot = mplot + delta;
    ctr = ctr + 1;
end