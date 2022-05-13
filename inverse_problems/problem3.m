clear all; 
close all; 
%% Preallocations
index = @(i) i + 1;                                                         % helper functions

n1    = 5;                                                                  % upper bounds for domain of discrete distributions of Xi
n2    = 10; 
n3    = 15;
N     = [n1 n2 n3];

y1    = 5.82;                                                               % measurements
y2    = 3.95;
y3    = 2.11;
Y     = [y1 y2 y3];

N1    = 0:n1;                                                               % domain of discrete distributions of Xi
N2    = 0:n2;
N3    = 0:n3;

sigma = sqrt(0.3);                                                          % stdev of the noise in the measurement
%% Posterior evaluation algortihm
post = zeros(length(N1), length(N2), length(N3));                           % posterior probability tensor

for k = N1                                                                  % compute values of the posterior distribution at every point
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
post = post/(sum(sum(sum(post))));                                          % normalization to unity
%% Conditional Mean Estimate
CM = zeros(length(N), 1);                                                   % conditional mean vector

for k = N1                                                                  % compute conditional mean of discrete posterior distribution
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

disp(CM)
%% MAP Estimate
[M,I]        = max(post,[],"all","linear");                                 % compute X_MAP as argmax p(X=[k l m)^T|Y=y)
[d1, d2, d3] = ind2sub(size(post),I);
map          = [d1-1, d2-1, d3-1];

disp(map)
%% Marginal posterior density of X1, X2
post12 = zeros(length(N1), length(N2));

for k = N1                                                                  % marginal distribution of X1, X2 is obtained by integrating out X3
    i = index(k);
    for l = N2
        j = index(l);
        for m = N3
            n = index(m);
            post12(i,j) = post12(i,j) + post(i,j,n);
        end
    end
end

visualize_marginal(post12, N1, N2, 1, "Marginal density of X1 and X2")      % visualize
visualize_marginal2d(post12, N1, N2, 2, "2D Marginal density of X1 and X2")
%% Marginal posterior density of X1, X3
post13 = zeros(length(N1), length(N3));

for k = N1                                                                  % marginal distribution of X1, X3 is obtained by integrating out X2
    i = index(k);
    for m = N3
        n = index(m);
        for l = N2
            j = index(l);
            post13(i,n) = post13(i,n) + post(i,j,n);
        end
    end
end

visualize_marginal(post13, N1, N3, 3, "Marginal density of X1 and X3")
visualize_marginal2d(post13, N1, N3, 4, "2D Marginal density of X1 and X3")
%% Marginal posterior density of X2, X3
post23 = zeros(length(N2), length(N3));

for l = N2                                                                  % marginal distribution of X2, X3 is obtained by integrating out X1
    j = index(l);
    for m = N3
        n = index(m);
        for k = N1
            i = index(k);
            post23(j,n) = post23(j,n) + post(i,j,n);
        end
    end
end

visualize_marginal(post23, N2, N3, 5, "Marginal density of X2 and X3")
visualize_marginal2d(post23, N2, N3, 6, "2D Marginal density of X2 and X3")

% On the basis of obtained results, it seems pretty likely that true value
% of x = [1 7 8]^T and difference with obtained MAP and CM estimates can be
% explained by presense of noise in the measurement values and ill posed
% nature of inverse reconstruction. 
%% Posterior function
function post = Posterior(Y, X, N, sigma)
    n          = length(X);
    
    amean      = sum(X)/n;
    gmean      = prod(X)^(1/n);
    hmean      = n/(1/X(1)+1/X(2)+1/X(3));

    likelihood = exp(-1/(2*sigma^2) * ...
                    (abs(Y(1) - amean)^2 + ...
                     abs(Y(2) - gmean)^2 + ...
                     abs(Y(3) - hmean)^2));

    prior      = 1/(2^sum(N))*prod(arrayfun(@(N, X) nchoosek(N,X), N, X));
    %prior      = 1/(2^sum(N))*nchoosek(N(1),X(1))*nchoosek(N(2),X(2))*nchoosek(N(3),X(3));
    
    post = likelihood * prior;
end
%% Visualization
function visualize_marginal(marg_dist, N, M, k, str)
    figure(k)
    stem3(N,M,marg_dist');
    title(str)
end

function visualize_marginal2d(marg_dist, N, M, k, str)
    figure(k)
    imagesc(marg_dist');
    set(gca,'YDir','normal');
    axis square
    title(sprintf(str), ...
            'FontSize', ...
            14);
end