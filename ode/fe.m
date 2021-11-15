clear all; 
close all; 

%(a) - Solving system using Forward Euler
x1 = @(x, alpha) -x^3/3 + x + alpha;
alpha1 = @(x, eps) -eps * x; 

%Initial conditions
eps = 0.001; 
%eps = 1;
x = 2; 
alpha = 2/3; 
%alpha = 4/3;

%Since we are using one of the most "inaccurate methods"
%lets choose small grid
h = 1/1000; 
t = 400; 

%storing all the approximations
x_approx = zeros(t, 1); 
alpha_approx = zeros(t, 1); 

for i = 1:t
    %record approximation of the step
    x_approx(i) = x; 
    alpha_approx(i) = alpha;
    %compute system values at the point
    fx = x1(x, alpha); 
    falpha = alpha1(x, eps);
    %estimate next argument
    x = x + h * fx; 
    %h * falpha
    alpha = alpha + h * falpha;
end


figure(1);
plot(1:t, x_approx); 
legend('x');
grid on; 
figure(2) 
plot(1:t, alpha_approx);
grid on; 
legend('alpha');

% Definitely we need to consider the restriction on step size. We can
% clearly observe that if we make a larger step size - the pattern of the
% function will be more "rough". And the function itself will contain
% more discontinuities. The choice of the step size mainly was motivated to
% capture the pattern of the reconstructed function. All in all, Euler
% method is not the most accurate and even with small step size it is badly
% estimating the function in question. Nevertheless, it is easy to
% implement. The choice of Euler method was made to demonstrate its ease of
% implementation and at the same time inaccuracy in estimation. 

% (c) - How sensitive is the system to small changes to eps and the initial
% values? 


%Lets change initial conditions and the interval
 
eps = 1;
x = 2; 
alpha = 4/3;


h = 1/1000; 
tt = 40000; 

%storing all the approximations
x_approx1 = zeros(tt, 1); 
alpha_approx1 = zeros(tt, 1); 

for i = 1:tt
    %record approximation of the step
    x_approx1(i) = x; 
    alpha_approx1(i) = alpha;
    %compute system values at the point
    fx = x1(x, alpha); 
    falpha = alpha1(x, eps);
    %estimate next argument
    x = x + h * fx; 
    %h * falpha
    alpha = alpha + h * falpha;
end


figure(3);
plot(1:tt, x_approx1); 
legend('x1');
grid on; 
figure(4) 
plot(1:tt, alpha_approx1);
grid on;
legend('alpha1'); 