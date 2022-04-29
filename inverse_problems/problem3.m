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
post = post/(sum(sum(sum(post))));                         % normalization
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
map          = [d1-1, d2-1, d3-1];
%% Computation of marginal densities 
post12 = zeros(length(N1), length(N2));

for k = N1
    i = index(k);
    for l = N2
        j = index(l);
        for m = N3
            n = index(m);
            post12(i,j) = post12(i,j) + post(i,j,n);
        end
    end
end

figure(1)
stem3(N1,N2,post12');
title("Marginal density of X1 and X2")

post13 = zeros(length(N1), length(N3));

for k = N1
    i = index(k);
    for m = N3
        n = index(m);
        for l = N2
            j = index(l);
            post13(i,n) = post13(i,n) + post(i,j,n);
        end
    end
end

figure(2)
stem3(N1,N3,post13');
title("Marginal density of X1 and X3")

post23 = zeros(length(N2), length(N3));

for l = N2
    j = index(l);
    for m = N3
        n = index(m);
        for k = N1
            i = index(k);
            post23(j,n) = post23(j,n) + post(i,j,n);
        end
    end
end

figure(3)
stem3(N2,N3,post23');
title("Marginal density of X2 and X3")
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
    %prior      = 1/(2^sum(N))*nchoosek(N(1),X(1))*nchoosek(N(2),X(2))*nchoosek(N(3),X(3));
    
    post = likelihood * prior;
end
%% Marginal density
function marg = Marginal(N1, N2, N3, post)
    marg = zeros(1);
    for k = N1
        i = index(k);
        for l = N2
            j = index(l);
            for m = N3
                n = index(m);
                post12(i,j) = post12(i,j) + post(i,j,n);
            end
        end
    end
end