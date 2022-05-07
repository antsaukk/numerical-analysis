clear all;
close all; 
%% Preallocations
index       = @(ind) ind-1;                                                 % setting helpers functions

U           = load("assignment2.mat");                                      % loading data from the file
Y           = U.Y;
w           = U.w;

N           = length(w);                                                    % setting the parameters
axis_length = sqrt(N);
h           = 1/axis_length;
x1          = linspace(0, 1-h, axis_length);
x2          = linspace(0, 1-h, axis_length);

%A           = Grid1(N, axis_length, x1, x2, Y, index);
A           = Grid2(N, axis_length, x1, x2, Y, index);

z           = zeros(size(w));                                               % estimate from alternating algorithm
gamma       = 0.1;                                                          % estimate for stdev

sigma       = 0.005;                                                        % standard deviation in of noise in the data
eps         = sigma * sqrt(N);                                              % Morozov criteria 

Aw = A'*w;                                                                  % precompute constant terms 
%% CGA

pcg             = zeros(N, 1);                                              % saving iterates of CGA
residuals_cg    = [];                                                       % saving residuals of CGA
k               = 1;                                                        % iterations counter

r               = Aw - A'*(A*pcg);                                          % first step of CGA
s               = r;
res             = norm(A*pcg - w);
residuals_cg(k) = res;

run             = res > eps;                                                % run criteria: until Morozov is satisfied

tic
while(run)
    AAs             = A'*(A*s);                                             % precompute iteration-constant terms to speed-up
        
                                                                            % Conjugate Gradient Algorithm
    alpha           = norm(r)^2/(s'*AAs);
    pcg             = pcg + alpha*s;
    rk              = r - alpha*AAs;
    beta            = (norm(rk)/norm(r))^2;
    s               = rk + beta*s;
    k               = k + 1;
    
    r               = rk;                                                   % change of variables: rk=r_k+1
        
    res             = norm(A*pcg - w);                                      % record residuals
    residuals_cg(k) = res;

    run             = res > eps;                                            % check run criteria
end
toc

visualize2d(pcg, axis_length, 1, k, 'Reconstruction by CGA with k %d');     % plot reconstruction
%% Alternating algorithm to find Z_MAP and gamma_MAP
kk          = 1;                                                            % iterations counter
gammas      = [];                                                           % recording residuals of gamma
gammas(1)   = gamma;                                                        % initial value of gamma

run = norm(A*z - w) > eps; %wrong                                           % stopping according to Morozov criteria                                      
norm(A*z - w)

tic
while(run) % run
    delta       = sigma^2/gamma^2;                                          % compute minimizer for Z_map using Tikhonov functional
    K           = [A; sqrt(delta) * eye(N)];
    v           = [w; zeros(N, 1)];
    z           = K\v;

    gamma       = norm(z)/sqrt(N);                                          % compute gamma from T'(z,gamma)=0
    run         = norm(A*z - w) > eps;
    
    kk          = kk + 1;
    gammas(kk)  = gamma;                                                    % record gamma
end
toc

figure(2)
hold on
plot(1:kk, gammas);
title('Prior stded as function of iterations')
legend('Prior stdev')
grid on
hold off

visualize2d(z, axis_length, 3, kk, 'Reconstruction via statistical inversion with k %d');
%% Function for generation Grid1
function A = Grid1(N, axis_length, x1, x2, Y, index)
    A = zeros(N, N);
    tic
    for k = 1:N
        for i = 1:axis_length
            for j = 1:axis_length
                x        = [x1(i) x2(j) 0];
                y        = Y(:,k);
                ind      = j+index(i)*axis_length;
                A(k,ind) = 1/norm(x - y);
            end
        end
    end
    A = 1/N * A;
    toc
end
%% Function for generation Grid2
function A = Grid2(N, axis_length, x1, x2, Y, index)
    tic
    %grid = [];
    grid = zeros(3, N);

    for i = 1:axis_length
        for j = 1:axis_length
            point                           = [1/80 + index(i)/axis_length
                                               1/80 + index(j)/axis_length 
                                               0];
            grid(:, j+index(i)*axis_length) = point;
            %grid = [grid, point];
        end
    end
    
    A = zeros(N, N);

    for k = 1:N
        for l = 1:N
            A(k,l) = 1/norm(grid(:,l) - Y(:, k));
        end
    end
    A = 1/N * A;

    toc
end
%% Visualization function
function visualize2d(pix, n, k, iter, str)
    Pix = reshape(pix, n, n);

    figure(k)
    imagesc([0,1], [0,1], Pix);
    set(gca,'YDir','normal');
    axis square
    title(sprintf('Reconstruction via statistical inversion with k %d', ...
            iter-1), ...
            'FontSize', ...
            14);
end