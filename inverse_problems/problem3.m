clear all; 
close all; 
%% Preallocations
index = @(i) i + 1;

n1    = 5; 
n2    = 10; 
n3    = 15;
N     = [n1 n2 n3];

y1    = 5.82;
y2    = 3.95;
y3    = 2.11;
Y     = [y1 y2 y3];

N1    = 0:n1; 
N2    = 0:n2;
N3    = 0:n3;

sigma = sqrt(0.3);
%% Posterior evaluation algortihm
post = zeros(length(N1), length(N2), length(N3));

for k = N1
    i = index(k);
    for l = N2
        j = index(l);
        for m = N3
            n           = index(m);
            X           = [k l m];
            post(i,j,n) = Posterior(Y, X, N, sigma);
        end
    end
end
%% Conditional Mean Estimate
CM = zeros(length(N), 1);

for k = N1
    i = index(k);
    for l = N2
        j = index(l);
        for m = N3
            n    = index(m);
            prob = post(i,j,n);
            CM   = CM + [k l m]'*prob;
        end
    end
end
%% MAP Estimate
[M,I]        = max(post,[],"all","linear");
[d1, d2, d3] = ind2sub(size(post),I);
map          = [d1, d2, d3];
%% Posterior function
function post = Posterior(Y, X, N, sigma)
    n          = length(X);
    
    amean      = sum(X)/n;
    gmean      = prod(X)^(1/n);
    hmean      = n/(1/X(1)+1/X(2)+1/X(3));

    likelihood = exp(-1/(2*sigma^2) * ...
                    ((Y(1) - amean)^2 + ...
                     (Y(2) - gmean)^2 + ...
                     (Y(3) - hmean)^2));

    prior      = 1/(2^sum(N))*prod(arrayfun(@(N, X) nchoosek(N,X), N, X));
                 %nchoosek(N(1),X(1)) * 
                 %nchoosek(N(2),X(2)) * 
                 %nchoosek(N(3),X(3));
    
    post = likelihood * prior;
end