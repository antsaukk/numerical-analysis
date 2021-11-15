clear all; 
%Taylor
%derivatives
xn1 = @(t, x) 1 + x^2 + t^3;
xn2 = @(t, x, x1) 2 * x * x1 + 3*t^2; 
xn3 = @(t, x, x1, x2) 2 * x1^2 + 2 * x * x2 + 6 * t; 
xn4 = @(x, x1, x2, x3) 6 * x1 * x2 + 2 * x * x3 + 6; 
%x + 1
x_next = @(x, x1, x2, x3, x4, h) x + h * x1 + h^2/2 * x2 + h^3/6 * x3 + h^4/24 * x4; 

%initial conditions
t = 0; 
x = 0;
h = 1/128; 
n = 128; 

%for debug
tt = zeros(128,1); 
ii = zeros(128,1);
xt = zeros(128,1);

%loop over x 
for i = 1:n
    xt(i) = x; 
    x1 = xn1(t, x); 
    x2 = xn2(t, x, x1); 
    x3 = xn3(t, x, x1, x2); 
    x4 = xn4(x, x1, x2, x3); 
    x = x_next(x, x1, x2, x3, x4, h); 
    t = i * h; 
end
 
taylor = x; 

%Euler

%initial conditions
t = 0; 
x = 0;
h = 1/128; 
n = 128;

xe = zeros(128,1); 

%loop over the x 
for i = 1:n
    xe(i) = x;
    x1 = xn1(t, x); 
    x = x + h*x1; 
    t = i * h; 
end

euler = x; 

figure(1);
plot(1:n, xt); 
hold on; 
plot(1:n, xe);
legend('taylor estimates', 'euler estimates');
grid on; 


