ebruary20clc;clear;clf; close all;
n  = 100;           %no of intervals
pn = 0; p0 = 0;     %boundaries
p = zeros([1 n]);   %alloc
p(1)=p0;
p(n)=pn;
pnew = p+1;
c = 1;
dx = 1/n;
yps = 10^-11; m=1; ctr = 0;
while m>yps 
    i = 1;
    pnew(1)=p(1);
    pnew(n)=p(n);
    %fprintf('%5.3f ',pnew(1));
    while i<n-1
        i = i + 1;
        pnew(i)=(dx*dx*f(i*dx)+p(i+1)+p(i-1))/2;
        %fprintf('%5.3f ',pnew(i));
    end
%     fprintf('%5.3f ',pnew(n));
%     fprintf('\n')
    m = max(abs(pnew-p));
%     if mod(ctr,100)==0;
%         plot(p); hold on; pause(0.05);
%     end
    p = pnew;
    ctr = ctr + 1;
end

x = 0:1/(n-1):1; y = exp(x)+(1-exp(1))*x-1; %analytic sol'n 
plot(x,y,'r'); hold on; 
plot(x,p,'k') %plots

i = 1;  maxi = 0; %finding max error 
while i < 101 
    if abs(p(i)-y(i)) > maxi
        maxi = abs(p(i)-y(i)); maxirel = maxi/abs(y(i)); index = i;
    end
    i  = i + 1;
end


function y = f(x)
     y = -1*exp(x);
end


% The Jacobi algorithm behaves appropriate, confidently being around the analytical solution. 
% The absolute error said to be 0.003300400349183 on [-0.2; 0] segment range which comes to be the 1.5615506114819% in the middle of the domain in solution.
% 
% Number of iterations, however turns out to be 27588 having epsilon 10^-10.
% Concluding  it is fair to declare that the Jacobi method is not the fastest algorithm to be used, especially in time limited situations.
% 
% Rauan Kelesbekov