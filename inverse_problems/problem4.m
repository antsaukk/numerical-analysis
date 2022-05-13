clear all;
close all;
%% Preallocations
T          = 0.02;                                                          % time instance
N          = 100;                                                           % number of discretization intervals
h          = pi/N;                                                          % interval width
x          = linspace(h, pi-h, N-1);                                        % grid discretization

B          = diag(ones(N-2, 1), -1)-2*eye(N-1)+diag(ones(N-2, 1), 1);       % approximation of 2nd order differential operator
B(N-1,N-1) = -1;
B          = 1/h^2 * B;

A          = expm(T*B);                                                     % matrix exponential operator to solve forward problem? Inverse crime?

Y          = load("assignment4.mat");                                       % loading data struct
y          = [Y.y; zeros(N-1, 1)];                                          % noisy measurement

sigma      = 0.5;                                                           % stdev of noise in measurement data
eps        = sigma * sqrt(N-1);                                             % acceptance level according to Morozov discrenpancy principle
stop       = 10^(-10);                                                      % stopping criteria

F_Tikhonov = zeros(1, length(Y.y));                                         % data structure to store Tikhonov functional solution

N_sam      = 5000;                                                           % number of samples
gamma      = 5;                                                             % value of gamma parameter
F_Gibbs    = zeros(1, length(Y.y));                                         % reconstruction from gibbs sampler
%% (a) Solution to heat equation via Minimization of Tikhonov functional

delta = 0.1;                                                                % starting value for delta
run   = true;
while(run)
    K     = [A; sqrt(delta)*eye(N-1)];                                      % K matrix for of Tikhonov functional
    xdh   = K\y;                                                            % current estimate for initial heat dist.
    f     = norm(A*xdh - y(1:N-1))^2 - eps^2;                               % discrenpancy function
    df    = 2*delta*xdh'*inv(A'*A+delta*eye(N-1))*xdh;                      % classical derivative of discrenpancy function
    delta = delta - f/df;                                                   % newton step
    run   = (abs(f) > stop);                                                % stopping criteria
end
K          = [A; sqrt(delta)*eye(N-1)];
F_Tikhonov = K\y;                                                           % reconstructed solution
%% (b) Solution to heat equation equation with truncated SVD and Morozov principle

[U L V] = svd(A);                                                           % svd
diagonal = diag(L);                                                         % singular values
F_TSVD = zeros(1, length(Y.y));                                             % data structure to store TSVD solution
residuals = [];                                                             % store residuals
y = Y.y; 

cutoff = 1;
for n = 1:length(diagonal)
    sv = [1./diagonal(1:n); zeros(N-n-1, 1)];                               % inverte singular values
    A_pinv_n = V * diag(sv) * U';                                           % compute truncated svd

    F_TSVD = A_pinv_n * y;                                                  % solve inverse problem 
    
    residual = norm(A*F_TSVD - y) - eps;                                    % compute Morozov Criteria
    residuals(n) = residual;
    if norm(A*F_TSVD - y) <= eps 
        break;
    end
    cutoff = cutoff + 1;
end
disp(cutoff)
%% Gibbs sampler and conditional mean
x_range   = linspace(0, 30, 100);
Zsamples  = zeros(N_sam, N-1);
Ztemp     = zeros(size(F_Gibbs));  

tic
for s = 1:N_sam
    for j = 1:(N-1)
        I_line = zeros(size(x_range));  %???                                     % take jth component and evaluate cond.dist. over line
        for i = 1:length(I_line)
            Ztemp     = F_Gibbs;
            Ztemp(j)  = x_range(i);
            I_line(i) = PosteriorIter(Ztemp, Y.y, A, gamma, sigma);
        end

        cdf              = cumsum(I_line);                                    % integrate condional density
        cdf              = cdf/cdf(end);                                         % normalization to 1
        tau              = rand;                                                 % sample from Unif(0,1)
        xi               = find(tau <= cdf, 1);                                  % inverse of cdf approximated numerically

        F_Gibbs(j)       = x_range(xi);
    end
    Zsamples(s, :) = F_Gibbs;
end
toc

postsamples_probs = zeros(1, N_sam);
for i = 1:N_sam
    postsamples_probs(i) = Posterior(Zsamples(i, :), Y.y, A, gamma, sigma);
end

normalized_post = postsamples_probs/sum(postsamples_probs);
ZCM = zeros(length(Y.y), 1);
for i = 1:N_sam
    ZCM = ZCM + Zsamples(i, :)'*normalized_post(i);
end

indexes_of_min_mix_axes = sample_histories(Zsamples);
indexes_of_min_mix_axes
%% Plots
figure(1);
hold on
plot(x, F_Tikhonov, 'LineWidth', 1);
title('Tikhonov reconstruction of initial heat distribution')
legend('Approximation F via minimization of TF')
grid on;
hold off

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

figure(5);
hold on
plot(x, ZCM, 'LineWidth', 1);
title('GS reconstruction of initial heat distribution')
legend('Approximation of F obtained by sampling')
grid on;
hold off

figure(6);
hold on
plot(x, A*F_Tikhonov, 'LineWidth', 1);
plot(x, A*F_TSVD, 'LineWidth', 1);
plot(x, A*ZCM, 'LineWidth', 1);
plot(x, Y.y, 'LineWidth', 1);
cap = sprintf('A*Tikhonov, A*TSVD, A*GS, Measurement');
title(cap, 'FontSize', 10)
legend('A*(Tikhonov solution)', ...
    'A*(TSVD solution)', ...
    'A*(GS Solution)', ...
    'Measurement')
grid on;
hold off

figure(7)
hold on
subplot(2,1,1);plot(Zsamples(:, indexes_of_min_mix_axes(1)));
title('Sample history smallest');
subplot(2,1,2);plot(Zsamples(:, indexes_of_min_mix_axes(2)));
title('Sample history largest');
hold off
%% Utilities
function post = Posterior(z, y, A, gamma, sigma)
    Ksi030 = @(z) z >= 0 && z <= 30;
    pi030  = prod(arrayfun(@(z) Ksi030(z), z));
    pr_j   = @(z, gamma) 1/(1 + gamma^2*z^2);
    gamma  = gamma*ones(size(z));
    prior  = prod(arrayfun(@(z, gamma) pr_j(z, gamma), z, gamma));
    
    likeli = @(y, A, z, sigma) exp(-1/(2*sigma^2)*norm(y-A*z')^2);
   
    post   = pi030 * prior * likeli(y, A, z, sigma);
end

function post = PosteriorIter(z, y, A, gamma, sigma)
    post   = 1;
    for i = 1:length(z)
        post = post * 1/(1 + gamma^2*z(i)^2);
    end
    post = post * exp(-1/(2*sigma^2)*norm(y-A*z')^2);
end

function post = Posterior1D(z, y, ai, gamma, sigma, index)
    ksi030 = @(z) z >= 0 && z <= 30;
    pi030  = prod(arrayfun(@(z) ksi030(z), z));

    zi     = z(index);
    yi     = y(index);

    prior  = 1/(1 + gamma^2*zi^2);

    likeli = exp(-1/(2*sigma^2)*abs(yi-ai*z')^2);
    
    post   = pi030 * prior * likeli;
end

function indexes = sample_histories(Z)
    ZM      = mean(Z);
    argmin  = find(ZM == min(ZM));
    argmax  = find(ZM == max(ZM));
    indexes = [argmin, argmax];
end