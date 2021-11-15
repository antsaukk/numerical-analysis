%BISECTION - idea given by Ovall p.49, reimplemented by @antsaukk 
f = @(x)(3*x-2)*(x*x+1);
a = 0; 
b = 1;
nmax=100; 
tol=1.0e-8;
fa=f(a); 
fb=f(b); 
if(fa*fb >= 0) 
    disp('Bad search interval'); 
    x=NaN; fx=NaN; n=0; d=NaN; 
    return; 
end
x = zeros(26, 1);
x(1)=(a+b)/2; 
fx=f(x(1)); 
n=0; 
d=(b-a)/2; 
i = 1;
while(n < nmax & abs(fx)>tol & d>tol) 
    if(fa*fx<0) 
        b=x(i);
        fb=fx; 
    else
        a=x(i);
        fa=fx; 
    end
    x(i+1)=(a+b)/2; 
    fx=f(x(i+1)); 
    n=n+1; 
    d=d/2; 
    i=i+1;
end
figure(2);
plot(0:n - 1,x(1:n),'ko-')
title('Convergence of the bisection method.')
legend('convergence of x_{n} approximations')
xlabel('iterations')
ylabel('x_{n}')