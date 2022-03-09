clear all;
close all;
%%
% (a) 
% Devise a numerical scheme for computing the heat distribution at the time
% T = 0.1 for parabolic heat equation.
% Derivation is provided on paper. Numerics:

% Discretization of grid
T = 0.1;
N = 100;
h = pi/N;
x = linspace(h, pi-h, N-1);

% Approximation of differential operator
G = diag(ones(N-2, 1), -1) - (2+h^2) * eye(N-1) + diag(ones(N-2, 1), 1);
G(1,1) = G(1,1) + 1;
G = 1/h^2 * G;

% matrix exponential operator to solve the forward problem
A = expm(T*G);

% (b) 
% Solve the inverse problem of reconstructing Fj, j = 1,2, 
% by using Tikhonov regularization via minimization of Tikhonov functional

% load measurements
U = load("ex2.mat");

% form the problem
sigma = 0.01;
eps = sigma * sqrt(N-1);
stop = 10^-15;

z1 = [U.y1; zeros(N-1, 1)];
z2 = [U.y2; zeros(N-1, 1)];
Z = [z1, z2];

% store answers
F = zeros(2, length(U.y1));
D = zeros(2, 1);

% Newton method
ind = 1;
for z = Z
    delta = 0.1;
    run = true;
    while(run)
        % K matrix for of Tikhonov functional
        K = [A; sqrt(delta)*eye(N-1)];

        % current estimate for initial heat dist.
        xdh = K\z; 
        
        % f(\delta)
        f = norm(A*xdh - z(1:N-1))^2 - eps^2;

        % f'(\delta)
        df = 2*delta*xdh'*inv(A'*A+delta*eye(N-1))*xdh;

        % newton step
        delta = delta - f/df;
        
        % stopping criteria
        run = (f > stop);
    end
    K = [A; sqrt(delta)*eye(N-1)];
    
    % reconstructed solution
    F(ind,:) = K\z;

    % delta chosen with Morozov discrenpancy principle
    D(ind,:) = delta;
    ind = ind + 1;
end

%%
% (c)
% One of the two measurement vectors corresponds to the initial heat 
% distribution 
% f(x) = 1/10 x^2(2x−π)(3x−2π)(x−π) 
% and the other to 
% f(x) = 1/10 x^2(2x−π)(3x−2π)(x−π)+2sin(8.5(x−π)). 
% Can you tell which one is which? Is the result surprising? 
% Give an intuitive explanation for this phenomenon?

f1 = 1/10*x.^2.*(2*x-pi).*(3*x-2*pi).*(x-pi);
f2 = f1 + 2*sin(8.5*(x-pi));

figure(1);
hold on
plot(x, F(1,:), 'LineWidth', 1);
plot(x, f1, 'LineWidth', 1);
title('Reconstruction of F1 init. heat distribution from y1')
legend('Approximation F1', 'Analytical f1')
grid on;
hold off

figure(2)
hold on
plot(x, F(2,:), 'LineWidth', 1);
plot(x, f2, 'LineWidth', 1);
title('Reconstruction of F2 init. heat distribution from y2')
legend('Approximation F2', 'Analytical f2')
grid on;
hold off 

figure(3)
hold on
plot(x, F(1,:), 'LineWidth', 1);
plot(x, F(2,:), 'LineWidth', 1);
title('Reconstruction of F1 and F2 init. heat distributions from y1 and y2')
legend('F1', 'F2')
grid on;
hold off 

figure(4)
hold on
plot(x, f1, 'LineWidth', 1);
plot(x, f2, 'LineWidth', 1);
title('Analytical form of initial heat distrbutions')
legend('f1', 'f2')
grid on;
hold off 