clear all;
close all;
%% Preallocations
dfdx         = @(q,x,y,i) -q*(y(i)-x(i))/norm(y-x')^2;                 % dfdx component of Jacobian
dfdq         = @(x,y) log(norm(y-x'));                                      % dfdq component of Jacobian
fi           = @(x1,x2,q1,q2,y) q1*log(norm(y-x1')) + q2*log(norm(y-x2'));  % function of the measurement

U            = load("assignment1.mat");                                     % unpack and set data and parameters
Y            = U.Y;
w            = U.w;
N            = length(w);
M            = 6;

x1           = [-0.25, 0];                                                  % initial locations and charges
x2           = [0.25, 0];                                                   % of the particle
q1           = 0;
q2           = 0;

% Since compoentns of noise vector are mutually independent,
% normally distributed with 0 mean and sigma stdev, we can
% estimate Morozov discrenpancy level as Expected value of the squared norm
% of noise vector. It is easy to check that with this choice we terminate
% the algorithm when norm of the residual is smaller than
% epsion=sigma*sqrt(N). The detailed derivation of Morozov epsilon-acceptance 
% level is presented in the solution of the problem 4.

sigma        = 0.05;                                                        % stdev of noise in the data
eps          = sigma * sqrt(N);                                             % Morozov level of discrenpancy

deltas       = [100 10 1 10^(-5)];                                          % experimental values of regularization parameter
%deltas       = [10^(-5)];
iter_bound   = 100000;

z            = zeros(N, 1);                                                 % preallocations of data structures required for algorithm
I            = eye(M);
lp           = zeros(M, 1);
Z            = zeros(length(deltas), M);

zdk          = [x1, x2, q1, q2]';                                           % pack initial values of interest in vector for Tikhonov algorithm
%% Levenberg-Marquardt Algorithm
k          = 1;                                                             % iteration counter

for delta = deltas                                              
    run  = true;                                                            % initial running criteria
    iter = 0;                                                               % evaluate convergence

    tic
    while(run && iter < iter_bound)                                         % either satisfy Morozov level or choose bound to terminate iteration manually
        % Linearization components
        Jfz  = Jacobian(x1, x2, q1, q2, Y, dfdx, dfdq);                     % compute Jacobian
        fzk  = F(x1, x2, q1, q2, Y, fi);                                    % compute vector valued f

        % components of Tikhonov functional
        K    = [Jfz; sqrt(delta)*I];                                        % form the operator                                     
        v    = [w - fzk; lp];                                               % the rhs of minimizing equation
        z    = K\v;                                                         % compute current iterate using backslash
        
        zdk  = zdk + z;                                                     % compute next iterate

        x1   = [zdk(1), zdk(2)];                                            % unpack values of locations 
        x2   = [zdk(3), zdk(4)];                                            % and charges of the particle
        q1   = zdk(5);                                                      % related to current iterate
        q2   = zdk(6);

        stop = norm(Jfz*z + fzk - w)^2;                                     % check the stopping criteria
        %stop
        run  = stop > eps^2;                                                % wrt Morozov principle
        iter = iter + 1;                                                    % increase iterations
    end
    toc
    
    visualize(zdk, Y, delta, k, iter);                                      % visualize 

    Z(k,:)        = zdk;                                                    % save solution obtained with current delta
    k             = k + 1;

    x1            = [-0.25, 0];                                             % reinit initial values of
    x2            = [0.25, 0];                                              % locations and charges of the 
    q1            = 0;                                                      % particle and pack in the vector
    q2            = 0;
    zdk           = [x1, x2, q1, q2]';
end

% From the obtained reconstruction plots for three choices of 
% delta = [100 10 1] we observe that the estimated location and charges are 
% only marginally different, meaning that for all three deltas the
% algorithm has converged to the same solution vector and marginal
% difference can be explained by the fact that for the different choice of
% delta obtained reconstruction enters set of Morozov discrenpancy
% level with different magnitude. What is different, however, is the
% convergence speed. From the reported computation times and iterations
% that algorithm took to converge, we observe that choice of delta=1 seems
% to be an optimal one from the convergence speed point of view. For the delta =
% 10^(-5) we run the algorithm and pipe the value of residual into the
% standard output. In this case, we observe that algorithm fails to attain
% Morozov discrenpancy level and get stuck with residual being around 6.14.
% In this case, we do the following trick: we stop the algorthm manually
% for some reasonably large iteration bound, say 100000 and plot the
% resulting reconstructed solution. As we can see from the plots, the
% result is numerical garbage meaning the for the delta = 10^(-5)
% algorithm diverges. It is not surprising, since numerically speaking
% 10^(-5) is small value and performing large number of matrix arithmetics
% likely results in highly amplified floating point errors, which in fact
% distort the reconstructed solution.
%% Visualization function
function visualize(z, Y, delta, k, iterations)
    x1   = [z(1); z(2); z(5)];                                              % unpack the values
    x2   = [z(3); z(4); z(6)];
    X    = [x1, x2];

    xl   = -1;                                                              % parameters for visualization of the boundary
    xu   =  1;
    yl   = -1;
    yu   =  1;
    x    = [xl, xu, xu, xl, xl];
    y    = [yl, yl, yu, yu, yl];
    
    figure(k)                                                               % plot generation
    plot(x, y, 'b-', 'LineWidth', 3);
    hold on
    for xpart = X
        plot(xpart(1), xpart(2), '-o', 'LineWidth', 3)
        text(xpart(1), xpart(2),sprintf('q=%.3f', xpart(3)))
        hold on
    end
    
    for i = 1:length(Y)
        plot(Y(1, i), Y(2, i), '-o', 'LineWidth', 3, 'Color', 'k')
        text(Y(1, i), Y(2, i), sprintf('  M-DEVICE %d', i))
        hold on
    end
    title(sprintf('Solution with delta %d and iterations to converge %d', ...
        delta, iterations), 'FontSize', 14);
    hold off
end
%% Jacobian function
function Jf = Jacobian(x1, x2, q1, q2, Y, dfdx, dfdq)
    Jf = zeros(8, 6);
    for i = 1:8
        jfi = [dfdx(q1,x1,Y(:,i), 1)                                        % compute each row of Jacobian
               dfdx(q1,x1,Y(:,i), 2)
               dfdx(q2,x2,Y(:,i), 1)
               dfdx(q2,x2,Y(:,i), 2)
               dfdq(x1,Y(:,i))
               dfdq(x2,Y(:,i))
               ];
        Jf(i,:) = jfi;
    end
end
%% Measurement function
function f = F(x1, x2, q1, q2, Y, fi)
    f = zeros(8, 1);
    for i = 1:8
        f(i) = fi(x1,x2,q1,q2,Y(:,i));
    end
end