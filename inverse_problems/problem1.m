clear all;
close all;
%% Preallocations
dfdx         = @(q,x,y,i) -q*(y(i)-x(i))/norm(y-x')^2;                      % dfdx component of Jacobian
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

sigma        = 0.05;                                                        % stdev of noise in the data
eps          = sigma * sqrt(N);                                             % Morozov level of discrenpancy

deltas       = [100 10 1];% 10^(-5)];                                       % experimental values of regularization parameter

z            = zeros(N, 1);                                                 % preallocations of data structures required for algorithm
I            = eye(M);
lp           = zeros(M, 1);
Z            = zeros(length(deltas), M);

zdk          = [x1, x2, q1, q2]';                                           % pack initial values of interest in vector for Tikhonov algorithm
%% Levenberg-Marquardt Algorithm
k = 1;                                                                      % iteration counter
for delta = deltas                                              
    run = true;                                                             % initial running criteria

    tic
    while(run)
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
        run  = stop > eps^2;                                                % wrt Morozov principle
    end
    toc

    Z(k,:) = zdk;                                                           % save solution obtained with current delta
    visualize(zdk, Y, k);                                                   % visualize 
    k = k + 1;

    x1           = [-0.25, 0];                                              % reinit initial values of
    x2           = [0.25, 0];                                               % locations and charges of the 
    q1           = 0;                                                       % particle and pack in the vector
    q2           = 0;
    zdk          = [x1, x2, q1, q2]';
end
toc
%% Visualization function
function visualize(z, Y, k)
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
    %plot(x1(1), x1(2), '-o', 'LineWidth', 3)
    %text(x1(1), x1(2),sprintf('q1=%.3f', q1))
    %hold on
    %plot(x2(1), x2(2), '-o', 'LineWidth', 3)
    %text(x2(1), x2(2),sprintf('q2=%.3f', q2))
    %hold on
    for i = 1:length(Y)
        plot(Y(1, i), Y(2, i), '-o', 'LineWidth', 3, 'Color', 'k')
        text(Y(1, i), Y(2, i),sprintf('  M-DEVICE %d', i))
        hold on
    end

    caption = sprintf('something k %d', k);
    title(caption, 'FontSize', 14);
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