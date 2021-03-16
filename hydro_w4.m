n=100;
curr=zeros(1,n)+1;
prev=zeros(1,n);
m=max(abs(curr-prev));
eps=0.000001;
C=1;
L=1;
dx = L/n;
dt=0.01;
alpha = prev;
beta = prev;
A=prev;
B=prev;
D=prev;
p1=1;
x = dx:dx:1;
every=10;
ctr =0;
while m>eps
    for i = 2:n-1
        curr(i)= -(C*dt/dx)*(prev(i)-prev(i-1))+prev(i);
    end
    curr(1)=1;
    curr(n)=curr(n-1);
    m=max(abs(curr-prev));
    prev=curr;
    if mod(ctr,every)==0
    plot(x,curr); hold on; pause(0.5)
    end
   ctr = ctr + 1;
end
