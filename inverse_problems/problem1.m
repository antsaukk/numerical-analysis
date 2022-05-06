clear all;
close all;
%% Preallocations
dfdx         = @(q,x,y,i) -q*(y(i)-x(i))/norm(y-x)^2; 
dfdq         = @(x,y) log(norm(y-x));
fi           = @(x1,x2,q1,q2,y) q1*log(norm(y-x1')) + q2*log(norm(y-x2'));

U            = load("assignment1.mat");
Y            = U.Y;
w            = U.w;
N            = length(w);
M            = 6;

x1           = [-0.25, 0];
x2           = [0.25, 0];
q1           = 0;
q2           = 0;

sigma        = 0.05;
eps          = sigma * sqrt(N);

z            = zeros(N, 1);
I            = eye(M);
lp           = zeros(M, 1);
deltas       = [100 10 1 10^(-5)];
Z            = zeros(length(deltas), M);
zdk          = [x1, x2, q1, q2]';
%% Levenberg-Marquardt Algorithm
k = 1;
for delta = deltas
    run = true;
    i = 1;
    tic
    while(run)
        % Linearization components
        Jfz  = Jacobian(x1, x2, q1, q2, Y, dfdx, dfdq);
        %Jfz
        fzk  = F(x1, x2, q1, q2, Y, fi);
        %fzk
        % Tikhonov function components
        %u    = Jfz*zdk' + w - fzk;
        K    = [Jfz; sqrt(delta)*I];
        %K
        %v    = [u; zdk'];
        %z    = K\v;
        v    = [w - fzk; lp];
        %v
        z    = K\v;
        %z
        zdk  = zdk + z;
        %zdk

        % new values unpacked
        x1   = [zdk(1), zdk(2)];
        x2   = [zdk(3), zdk(4)];
        q1   = zdk(5);
        q2   = zdk(6);

        % stopping criteria
        %stop = fzk + Jfz*(z - zdk') - w; %Matrix? 
        %run  = norm(stop) > eps;

        stop = norm(Jfz*z + fzk - w)^2;
        %stop
        run  = stop > eps^2;
        %i = i + 1;
        %if i == 4
        %    break;
        %end
    end
    toc
    Z(k,:) = zdk;
    k = k + 1;

    x1           = [-0.25, 0];
    x2           = [0.25, 0];
    q1           = 0;
    q2           = 0;
    zdk          = [x1, x2, q1, q2]';
end
toc
%% Visualizations

%% Jacobian function
function Jf = Jacobian(x1, x2, q1, q2, Y, dfdx, dfdq)
    Jf = zeros(8, 6);
    for i = 1:8
        %jfi = [dfdx(q1,x1,Y(:,i), 1) %@(q,x,y,i) -q1*(Y(1,i)-x1(1))/norm(Y(:,i)-x1)^2
        %       dfdx(q1,x1,Y(:,i), 2)
        %       dfdx(q2,x2,Y(:,i), 1)
        %       dfdx(q2,x2,Y(:,i), 2)
        %       dfdq(x1,Y(:,i))
        %       dfdq(x2,Y(:,i))
        %       ];
        jfi = [-q1*(Y(1,i)-x1(1))/norm(Y(:,i)-x1')^2
               -q1*(Y(2,i)-x1(2))/norm(Y(:,i)-x1')^2
               -q2*(Y(1,i)-x2(1))/norm(Y(:,i)-x2')^2
               -q2*(Y(2,i)-x2(2))/norm(Y(:,i)-x2')^2
               log(norm(Y(:,i)-x1'))
               log(norm(Y(:,i)-x2'))
               ];
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