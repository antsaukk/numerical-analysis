clear all;
close all; 
%% Preallocations
index       = @(ind) ind-1;                                                 % define helpers functions

U           = load("assignment2.mat");                                      % load data from the file
Y           = U.Y;
w           = U.w;

N           = length(w);                                                    % define the parameters
axis_length = sqrt(N);
h           = 1/axis_length;
x1          = linspace(0, 1-h, axis_length);
x2          = linspace(0, 1-h, axis_length);

%A           = Grid1(N, axis_length, x1, x2, Y, index);
A           = Grid2(N, axis_length, Y, index);                              % compute approximation of integral operator

z           = zeros(size(w));                                               % initial estimate for alternating algorithm
gamma       = 0.1;                                                          % initial estimate for stdev

% Since compoentns of noise vector are mutually independent,
% normally distributed with 0 mean and sigma stdev, we can
% estimate Morozov discrenpancy level as Expected value of the squared norm
% of noise vector. It is easy to check that with this choice we terminate
% the algorithm when norm of the residual is smaller than
% epsion=sigma*sqrt(N). The detailed derivation of Morozov epsilon-acceptance 
% level is presented in the solution of the problem 4.

sigma       = 0.005;                                                        % standard deviation of noise in the data
eps         = sigma * sqrt(N);                                              % Morozov criteria


Aw = A'*w;                                                                  % precompute constant terms 
%% CGA

pcg             = zeros(N, 1);                                              % save iterates of CGA
residuals_cg    = [];                                                       % save residuals of CGA
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
    r               = rk;                                                  
        
    res             = norm(A*pcg - w);                                      % record residuals
    residuals_cg(k) = res;

    run             = res > eps;                                            % check stopping criteria
end
toc

visualize2d(pcg, axis_length, 1, k, 'Reconstruction by CGA with k %d');     % plot reconstruction obtained by Conjugate Gradient Algorithm
%% Alternating algorithm to find Z_MAP and gamma_MAP
kk          = 1;                                                            % iterations counter
gammas      = [];                                                           % recording residuals of gamma
gammas(1)   = gamma;                                                        % initial value of gamma

% Since iterate of Z_MAP and gamma_MAP are dependent, it implies that it is
% enough to look only the at convergence Z_MAP, since as Z_MAP will be in
% the ballpark of Morozov criteria, the algorithm will stop updating it to
% avoind overfitting and it implies that gamma wont be updated either as it
% is derived as function of Z_MAP.
run = norm(A*z - w) > eps;                                                  % stopping according to Morozov criteria                                      
norm(A*z - w)

tic
while(run) % how to stop?
    delta       = sigma^2/gamma^2;                                          % fix gamma and compute minimizer for Z_map 
    K           = [A; sqrt(delta) * eye(N)];                                % as least square solution of derived equation
    v           = [w; zeros(N, 1)];
    z           = K\v;

    gamma       = norm(z)/sqrt(N);                                          % compute gamma from T'(z,gamma)=0 with value of current iterate for Z_map 
    run         = norm(A*z - w) > eps;                                      % check the Morozov criteria
    
    kk          = kk + 1;
    gammas(kk)  = gamma;                                                    % record current value of gamma iterate
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
function A = Grid2(N, axis_length, Y, index)
    tic
    grid        = zeros(3, N);
    inner_point = 1/(axis_length*2);
    %disp(inner_point)
    for i = 1:axis_length
        for j = 1:axis_length
            point                           = [inner_point + index(i)/axis_length
                                               inner_point + index(j)/axis_length 
                                               0];
            grid(:, j+index(i)*axis_length) = point;
        end
    end
    
    A = zeros(N, N);

    for k = 1:N
        for l = 1:N
            A(k,l) = 1/norm(grid(:, l) - Y(:, k));
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
    title(sprintf(str, ...
            iter-1), ...
            'FontSize', ...
            14);
end