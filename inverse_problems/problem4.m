clear all;
close all;
%% Preallocations
T = 0.02;                                                                   % time instance
N = 100;                                                                    % number of discretization intervals
h = pi/N;                                                                   % interval width
x = linspace(h, pi-h, N-1);                                                 % grid discretization

B = diag(ones(N-2, 1), -1) - 2 * eye(N-1) + diag(ones(N-2, 1), 1);          % approximation of 2nd order differential operator
B(N-1,N-1) = -1;
B = 1/h^2 * B;

A = expm(T*B);                                                              % matrix exponential operator to solve forward problem? Inverse crime?

Y = load("assignment4.mat");                                                % loading data struct
y = [Y.y; zeros(N-1, 1)];                                                   % noisy measurement

sigma = 0.5;                                                                % stdev of noise in measurement data
eps = sigma * sqrt(N-1);                                                    % acceptance level according to Morozov discrenpancy principle
stop = 10^(-10);                                                            % stopping criteria

F_Tikhonov = zeros(1, length(Y.y));                                         % to store estimated initial heat distribution
%% (a) Solution to heat equation via Minimization of Tikhonov functional

delta = 0.1;                                                                % starting value for delta
run = true;
while(run)
    K = [A; sqrt(delta)*eye(N-1)];                                          % K matrix for of Tikhonov functional
    xdh = K\y;                                                              % current estimate for initial heat dist.
    f = norm(A*xdh - y(1:N-1))^2 - eps^2;                                   % discrenpancy function
    df = 2*delta*xdh'*inv(A'*A+delta*eye(N-1))*xdh;                         % classical derivative of discrenpancy function
    delta = delta - f/df;                                                   % newton step
    run = (abs(f) > stop);                                                  % stopping criteria
end
K = [A; sqrt(delta)*eye(N-1)];
F_Tikhonov = K\y;                                                           % reconstructed solution

figure(1);
hold on
plot(x, F_Tikhonov, 'LineWidth', 1);
title('Tikhonov reconstruction of initial heat distribution')
legend('Approximation F via minimization of TF')
grid on;
hold off
%% Solution to heat equation equation with truncated SVD and Morozov principle

[U L V] = svd(A);                                                           % svd
diagonal = diag(L);                                                         % singular values
F_TSVD = zeros(1, length(Y.y));
residuals = [];
y = Y.y; 

for n = 1:length(diagonal)
    sv = [1./diagonal(1:n); zeros(N-n-1, 1)];                               % inverting singular values
    A_pinv_n = V * diag(sv) * U';                                           % truncated svd

    F_TSVD = A_pinv_n * y;                                                  % solve inverse problem
    residual = norm(A*F_TSVD - y) - eps;                                    % compute Morozov Criteria
    residuals(n) = residual;
    if norm(A*F_TSVD - y) <= eps 
        break;
    end
end
figure(2);
hold on
plot(x, F_TSVD, 'LineWidth', 1);
title('TSVD reconstruction of initial heat distribution')
legend('Approximation of F w. TSVD')
grid on;
hold off

figure(3);
hold on
plot(x, F_Tikhonov, 'LineWidth', 1);
plot(x, F_TSVD, 'LineWidth', 1);
cap = sprintf('Tikhonov vs TSVD. Error: %.3f', norm(F_TSVD-F_Tikhonov));
title(cap, 'FontSize', 10)
legend('Tikhonov solution', 'TSVD solution')
grid on;
hold off

figure(4);
hold on
plot(abs(residuals), 'LineWidth', 1);
title('Residual ')
legend('||Ax - y||')
grid on;
hold off