clear all;
close all;
%% Preallocations
dfdx         = @(q,x,y,i) -q*(y(i)-x(i))/norm(y-x)^2; 
dfdq         = @(x,y) log(norm(y-x));
fi           = @(x1,x2,q1,q2,y) q1*log(norm(y-x1)) + q2*log(norm(y-x2));

U            = load("assignment1.mat");
Y            = U.Y;
w            = U.w;
N            = length(w);
M            = 6;

x1           = [-0.25 0];
x2           = [0.25 0];
q1           = 0;
q2           = 0;

sigma        = 0.05;
eps          = sigma * sqrt(N);

z            = zeros(N, 1);
I            = eye(M);
deltas       = [100 10 1 10^(-5)];
Z            = zeros(length(deltas), M);
%% Levenberg-Marquardt Algorithm
k = 1;
for delta = deltas
    run = true;
    while(run)
        % Linearization components
        Jfz  = Jacobian(x1, x2, q1, q2, Y, dfdx, dfdq);
        zdk  = [x1 x2 q1 q2];
        fzk  = F(x1, x2, q1, q2, Y, fi);
        
        % Tikhonov function components
        u    = Jfz*zdk' + w - fzk;
        K    = [Jfz; sqrt(delta)*I];
        v    = [u; zdk'];
        z    = K\v;
        
        % new values unpacked
        x1   = [z(1) z(2)];
        x2   = [z(3) z(4)];
        q1   = z(5);
        q2   = z(6);
        % stopping criteria
        stop = fzk + Jfz*(z - zdk') - w;
        run  = abs(stop) > eps;
    end
    Z(k,:) = z;
    k = k + 1;
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