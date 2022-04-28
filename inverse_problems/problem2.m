clear all;
close all; 
%% Preallocations
N = 1600; 
h = 1/N;
x1 = linspace(0, 1-h, N);
x2 = linspace(0, 1-h, N);
%[X, Y] = meshgrid(x1, x2);
U = load("assignment2.mat"); 
Y = U.Y;
w = U.w;

A = zeros(N, N);
tic
for i = 1:N
    for j = 1:N
        x = [x1(i) x2(j) 0];
        y = Y(:,i);
        A(i,j) = 1/norm(x - y);
    end
end
toc
A = 1/N^2 * A;

sigma = 0.005; % standard deviation in of noise in the data
eps = sqrt(sigma^2 * N^2); % Morozov criteria 

Aw = A'*w; % precompute constant terms 
%% CGA

pcg = zeros(N, 1);
residuals_cg = [];
k = 1; 

% first step
r = Aw - A'*(A*pcg);
s = r;
res = norm(A*pcg - w);
residuals_cg(k) = res;

run = res > eps; % run criteria: until Morozov is satisfied

tic
while(run)
    AAs = A'*(A*s); % precompute to speed-up
        
    % Conjugate Gradient Algorithm
    alpha = norm(r)^2/(s'*AAs);
    pcg = pcg + alpha*s;
    rk = r - alpha*AAs;
    beta = (norm(rk)/norm(r))^2;
    s = rk + beta*s;
    k = k + 1;
    
    r = rk; % change of variables: rk=r_k+1
        
    res = norm(A*pcg - w); % record residuals
    residuals_cg(k) = res;
        
    run = res > eps; % run criteria
end
toc

Pcg = reshape(pcg, 40, 40);
figure(1)
imagesc(Pcg);
set(gca,'YDir','normal');
axis square
caption = sprintf('Reconstruction with k %d', k-1);
title(caption, 'FontSize', 14);