clear all;
close all;
%% Preallocations
dfdx         = @(q,x,y,i) -q*(y(i)-x(i))/norm(y-x')^2;                      % dfdx component of Jacobian
dfdq         = @(x,y) log(norm(y-x'));                                      % dfdq component of Jacobian
fi           = @(x1,x2,q1,q2,y) q1*log(norm(y-x1')) + q2*log(norm(y-x2'));  % function

U            = load("assignment1.mat");                                     % unpack data and parameters
Y            = U.Y;
w            = U.w;
N            = length(w);
M            = 6;

x1           = [-0.25, 0];                                                  % initial locations and charges
x2           = [0.25, 0];
q1           = 0;
q2           = 0;

sigma        = 0.05;                                                        % stdev of noise in the data
eps          = sigma * sqrt(N);                                             % Morozov level of discrenpancy

z            = zeros(N, 1);                                                 % preallocations for algorithm
I            = eye(M);
lp           = zeros(M, 1);
deltas       = [100 10 1];% 10^(-5)];                                       % set of values for delta parameter
Z            = zeros(length(deltas), M);

zdk          = [x1, x2, q1, q2]';                                           % initial values packed in vector
%% Levenberg-Marquardt Algorithm
k = 1;
for delta = deltas
    run = true;
    tic
    while(run)
        % Linearization components
        Jfz  = Jacobian(x1, x2, q1, q2, Y, dfdx, dfdq);                     % compute Jacobian
        fzk  = F(x1, x2, q1, q2, Y, fi);                                    % compute vector valued f

        % Tikhonov function components
        K    = [Jfz; sqrt(delta)*I];                                        
        v    = [w - fzk; lp];
        z    = K\v;                                                         % compute current iterate
        
        zdk  = zdk + z;                                                     % compute next iterate

        % new values unpacked
        x1   = [zdk(1), zdk(2)];
        x2   = [zdk(3), zdk(4)];
        q1   = zdk(5);
        q2   = zdk(6);

        % stopping criteria
        stop = norm(Jfz*z + fzk - w)^2;
        run  = stop > eps^2;
    end
    toc

    Z(k,:) = zdk;                                                           % solution with current delta
    visualize(zdk, k);
    k = k + 1;

    x1           = [-0.25, 0];                                              % reinit initial values
    x2           = [0.25, 0];
    q1           = 0;
    q2           = 0;
    zdk          = [x1, x2, q1, q2]';
end
toc
%% Visualization function
function visualize(z, k)
    x1   = [z(1), z(2)];                                                    % unpack the values
    x2   = [z(3), z(4)];
    q1   = z(5);
    q2   = z(6);

    xl   = -1;                                                              % parameters for boundary visualization
    xu   =  1;
    yl   = -1;
    yu   =  1;
    x    = [xl, xu, xu, xl, xl];
    y    = [yl, yl, yu, yu, yl];
    
    figure(k)                                                               % plot generation
    plot(x, y, 'b-', 'LineWidth', 3);
    hold on
    plot(x1(1), x1(2), '-o', 'LineWidth', 3)
    text(x1(1), x1(2),sprintf('q1=%.3f', q1))
    hold on
    plot(x2(1), x2(2), '-o', 'LineWidth', 3)
    text(x2(1), x2(2),sprintf('q2=%.3f', q2))
    caption = sprintf('something k %d', k);
    title(caption, 'FontSize', 14);
    hold off
end
%% Jacobian function
function Jf = Jacobian(x1, x2, q1, q2, Y, dfdx, dfdq)
    Jf = zeros(8, 6);
    for i = 1:8
        jfi = [dfdx(q1,x1,Y(:,i), 1)
               dfdx(q1,x1,Y(:,i), 2)
               dfdx(q2,x2,Y(:,i), 1)
               dfdx(q2,x2,Y(:,i), 2)
               dfdq(x1,Y(:,i))
               dfdq(x2,Y(:,i))
               ];
        %jfi = [-q1*(Y(1,i)-x1(1))/norm(Y(:,i)-x1')^2
        %       -q1*(Y(2,i)-x1(2))/norm(Y(:,i)-x1')^2
        %       -q2*(Y(1,i)-x2(1))/norm(Y(:,i)-x2')^2
        %       -q2*(Y(2,i)-x2(2))/norm(Y(:,i)-x2')^2
        %       log(norm(Y(:,i)-x1'))
        %       log(norm(Y(:,i)-x2'))
        %       ];
        Jf(i,:) = jfi;
    end
end
%% Operator function
function f = F(x1, x2, q1, q2, Y, fi)
    f = zeros(8, 1);
    for i = 1:8
        f(i) = fi(x1,x2,q1,q2,Y(:,i));
    end
end