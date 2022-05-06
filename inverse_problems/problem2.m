clear all;
close all; 
%% Preallocations
index       = @(ind) ind-1;

N           = 1600;
axis_length = sqrt(N);
h           = 1/axis_length;
x1          = linspace(0, 1-h, axis_length);
x2          = linspace(0, 1-h, axis_length);
U           = load("assignment2.mat"); 
Y           = U.Y;
w           = U.w;

%A           = Grid1(N, axis_length, x1, x2, Y, index);
A           = Grid2(N, axis_length, x1, x2, Y, index);

z = zeros(size(w));                                                         % estimate from alternating algorithm
gamma = 0.1;                                                                % estimate for stdev

sigma = 0.005;                                                              % standard deviation in of noise in the data
eps = sigma * sqrt(N);                                                      % Morozov criteria 

Aw = A'*w;                                                                  % precompute constant terms 
%% CGA

pcg = zeros(N, 1);
residuals_cg = [];
k = 1; 

% first step
r = Aw - A'*(A*pcg);
s = r;
res = norm(A*pcg - w);
residuals_cg(k) = res;

run = res > eps;                                                            % run criteria: until Morozov is satisfied

tic
while(run)
    AAs             = A'*(A*s);                                             % precompute to speed-up
        
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

    run = res > eps; % run criteria
end
toc

Pcg = reshape(pcg, 40, 40);
figure(1)
imagesc([0,1], [0,1], Pcg);
set(gca,'YDir','normal');
axis square
caption = sprintf('Reconstruction with k %d', k-1);
title(caption, 'FontSize', 14);
%% Alternating algorithm to find MAP
kk = 1;
gammas = [];
gammas(1) = gamma;

run = norm(A*z - w) > eps && abs(sigma-gamma) > 0.001; %wrong
norm(A*z - w)
tic
while(kk < 2) % run
    % find z
    delta       = sigma^2/gamma^2;
    K           = [A; sqrt(delta) * eye(N)];
    v           = [w; zeros(N, 1)];
    z           = K\v;

    % find gamma
    gamma       = sqrt(2/N) * sigma * norm(z);
    run         = norm(A*z - w) > eps; %wrong
    norm(A*z-w)^2
    kk          = kk + 1;
    gammas(kk)  = gamma;
end
toc

figure(2)
hold on
plot(1:kk, gammas);
title('Prior stded as function of iterations')
legend('Prior stdev')
grid on
hold off

Z = reshape(z, 40, 40);
figure(3)
imagesc([0,1], [0,1], Z);
set(gca,'YDir','normal');
axis square
caption = sprintf('Reconstruction via statistical inversion with k %d', kk-1);
title(caption, 'FontSize', 14);
%% Function for generation Grid1
function A = Grid1(N, axis_length, x1, x2, Y, index)
    A = zeros(N, N);
    tic
    for k = 1:N
        for i = 1:axis_length
            for j = 1:axis_length
                x   = [x1(i) x2(j) 0];
                y   = Y(:,k);
                ind = j+index(i)*axis_length;
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
            %size(point)
            %j + index(i)*axis_length
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