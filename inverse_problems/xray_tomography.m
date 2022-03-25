clear all; 
close all; 
%% Test
U = load("ex3.mat");

% code to produce sinogram
x = [zeros(1,20) ones(1,40) zeros(1,40)];
y = [zeros(25,1); ones(40,1); zeros(35,1)];
B = y*x;
figure(1)
imagesc(B), axis square, colormap gray
caption = sprintf('Test Image 1');
title(caption, 'FontSize', 14);

z = U.A*B(:);
Z = reshape(z,70,70); % corresponding sinogram
figure(2)
imagesc(Z), axis square, colormap gray
caption = sprintf('Test Image 2');
title(caption, 'FontSize', 14);
%% Preallocations

% Tomography matrix and Sinogram
A = U.A;
y = U.S(:);

% precompute constant terms 
Ay = A'*y;

% AA = A'*A; % precompute matrix x matrix product

% noise and Morozov criteria
sigma = 0.015;
eps = sqrt(sigma^2 * 70^2);

iter = 3;
%% (a) 
% Reconstruct the object corresponding to the sinogram S by solving Ax = y 
% with Landweber-Fridman iteration method
for beta = [3, 10, 0.001]
    % init vector for solution
    xlw = zeros(U.N^2, 1);

    % record the residuals
    residuals_lw = [];
    k = 1;
    res = norm(y-A*xlw);
    residuals_lw(k) = res;

    % stopping criteria
    run = res > eps;
    while(run)
        k = k + 1;
    
        % new Landweber-Fridman iterate
        xlw = xlw + beta * (Ay - A'*(A*xlw));
    
        % residual
        res = norm(y-A*xlw);
        residuals_lw(k) = res;
       
        % stopping criteria
        run = res > eps;
        %res
    end
    % visualize reconstruction 
    Xlw = reshape(xlw,U.N,U.N);
    figure(iter)
    imagesc(Xlw), axis square, colormap gray
    caption = sprintf('Reconstruction with beta %d', beta);
    title(caption, 'FontSize', 14);

    iter = iter + 1;

    % visualize residual as function of iterations
    figure(iter);
    hold on
    plot(0:(k-1), residuals_lw, 'LineWidth', 1);
    caption_fig = sprintf('Residual of Landweber-Fridman vs k=%d with beta=%d', k-1, beta);
    title(caption_fig)
    legend('Residual of LW')
    grid on;
    hold off

    iter = iter + 1;
end

% Why beta = 10 is bad choice? Using Theorem from lecture notes we check
% that 0 < beta < 2/lambda_1^2, where lambda_1 is largest singular value
L = svds(A);
fprintf('Upper bound for beta to satisfy Theorem %d \n', 2/(max(L)^2));

% So it is clear that condition of the Theorem is not satified for beta =
% 10 and satisfied for beta = 3 and beta = 0.001, which is also confirmed
% by the plots with the reconstructed images: they look reasonable for
% beta=3 and beta=0.001, but useless for beta=10.
%% (b) 
% Solve reconstruction problem by applying conjugate gradient method to 
% normal equation A^TAx = A^Ty with x0 = 0.
N = 1000;
for index = [1, 2]
    % preallocation of solution, residuals, iteration variable
    xcg = zeros(U.N^2, 1);
    residuals_cg = [];
    k = 1;
    
    % first step
    r = Ay - A'*(A*xcg);
    s = r;
    res = norm(A*xcg - y);
    residuals_cg(k) = res;

    % run criteria:
    % until Morozov is satisfied
    % or for 1000 iterations
    if index == 1
        run = res > eps;
    else 
        run = k <= N;
    end

    while(run)
        % precompute to speed-up
        AAs = A'*(A*s);
        
        % Conjugate Gradient Algorithm
        alpha = norm(r)^2/(s'*AAs);
        xcg = xcg + alpha*s;
        rk = r - alpha*AAs;
        beta = (norm(rk)/norm(r))^2;
        s = rk + beta*s;
        k = k + 1;
        
        % change of variables: rk=r_k+1
        r = rk;
        
        % record residuals
        res = norm(A*xcg - y);
        residuals_cg(k) = res;
        
        % run criteria
        if index == 1
            run = res > eps;
        else 
            run = k <= N;
        end
    end
    
    % visualize reconstruction 
    Xcg = reshape(xcg,U.N,U.N);
    figure(iter)
    imagesc(Xcg), axis square, colormap gray
    caption = sprintf('Reconstruction with k %d', k-1);
    title(caption, 'FontSize', 14);

    iter = iter + 1;
    
    % residual error as function of k
    figure(iter);
    hold on
    plot(0:(k-1), residuals_cg, 'LineWidth', 1);
    caption_fig_cg = sprintf('Residual of Conjugate Gradient Method vs k=%d', k-1);
    title(caption_fig_cg)
    legend('Residual of CGM')
    grid on;
    hold off

    iter = iter + 1;
end

% As observable, if Conjugate Gradient Algorithm is run without early 
% stopping, poor reconstruction is obtained since the algorithm starts to 
% fit and reconstruct the numerical noise into the image.