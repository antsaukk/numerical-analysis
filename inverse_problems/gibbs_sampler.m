clear all;
close all;
% Gibbs sampler to reconstruct posterior probability distribution of the
% particle moving in the unit disk
%% Preallocations
f        = @(x,p) 1./abs(x-p);                                                   % charge potential function
post     = @(x,q,v,s,p) (abs(x)<=1) * exp(-norm(v-q*f(x,p))^2/(2*s^2));          % sampling posterior distribution

N        = 5000;                                                                 % number of samples
dim      = 2;                                                                    % number of dimensions

q        = 1;                                                                    % charge of particle
x0       = 0.5 + 0.4i;                                                           % location of the particle
sigma    = 0.15;                                                                 % measurement variance
n        = 3;                                                                    % number of sensors
theta    = linspace(0, 2*pi, n+1);                                               % uniform intervals
p        = exp(1i*theta(1:end-1));                                               % sensor locations
sigma_n  = sigma*max(abs(f(x0,p)));                                              % noise in measurement
v        = q*f(x0,p) + sigma_n*randn(n, 1);                                      % measurement with noise

t        = 0:pi/50:2*pi;                                                         % points to generate boundary of unit disk

xgp      = -1:0.001:1;                                                           % x-grid points
ygp      = xgp;                                                                  % y-grid points
[X, Y]   = meshgrid(xgp, ygp);                                                   % grid
compl_pl = complex(X, Y);                                                        % discretized complex plane
XYgp_com = [complex(xgp,0); complex(0,ygp)];                                     % complex absciss and ordinate for samplig from cond. density
%% Gibbs Sampler

sample_histories = zeros(N, 2);                                                  % record sample histories
posterior        = zeros(size(compl_pl));                                        % value of posterior at every grid point
C_density        = zeros(size(XYgp_com(1,:)));                                   % values to evalue conditional density over integration line

xk = 0 + 0i;
tic
for k = 1:N
    %y = 0 + 0i;
    for j = [1, 2]                                                               
        if j == 1                                                                % take jth component of xk-sampled point and add up to required axis to form integration line
            I_line = XYgp_com(j, :) + imag(xk)*1i;
        else
            I_line = XYgp_com(j, :) + real(xk);
        end
        
        for i = 1:length(I_line)                                                 % evaluate conditional density over integration line
            C_density(i) = post(I_line(i), q, v, sigma_n, p);
        end
        %C_density = arrayfun(@(I_line) post(I_line, q, v, sigma_n, p), I_line);

        cdf     = cumsum(C_density);                                             % integrate condional density
        cdf     = cdf/cdf(end);                                                  % normalization to 1
        tau     = rand;                                                          % sample from Unif(0,1)
        xi      = find(tau <= cdf, 1);                                           % inverse of cdf approximated numerically
        
        xk       = I_line(xi);
    end
    
    sample_histories(k, 1) = real(xk);
    sample_histories(k, 2) = imag(xk);
end
toc
%% Posterior density

for i = 1:length(Y)                                                              % evaluate posterior                                             
    for j = 1:length(X)
        posterior(i,j) = post(compl_pl(i,j),q,v,sigma_n,p); 
    end
end
%% Autocovariances

ac_x = autocovariance(sample_histories(:,1));
ac_y = autocovariance(sample_histories(:,2));
%% Plots
figure(1)
imagesc([-1,1], [-1,1], posterior);
set(gca,'YDir','normal');
axis image
hold on
plot(cos(t), sin(t), 'k');
plot(p, '*');
plot(x0, '+');
title('Posterior');
hold off

figure(2)
hold on
scatter(sample_histories(:,1), sample_histories(:,2),25,"red");
plot(cos(t), sin(t), 'k');
plot(p, '*', 'Color', 'b');
plot(x0, '+', 'Color', 'k');
cap = sprintf('Sample scatterplot');
title(cap, 'FontSize', 10);
hold off

figure(3)
hold on
subplot(2,1,1);plot(sample_histories(:,1));
title('Sample history x');
subplot(2,1,2);plot(sample_histories(:,2));
title('Sample history y');
hold off

L = 100;
figure(4)
hold on
subplot(2,1,1);plot(0:L, ac_x(1:(L+1)));
legend('ac of x')
title('Autocovariance of component x generated by Gibbs sampler');
subplot(2,1,2);plot(0:L, ac_y(1:(L+1)));
legend('ac of y')
title('Autocovariance of component y generated by Gibbs sampler');
hold off
%% Autocovariance function
function ac = autocovariance(x)
    mn       = mean(x);
    centered = x - mn;
    gamma0   = mean(centered.^2);
    N        = length(x);
    K        = 1:(N-1);
    ac       = zeros(1, N-1);
    for k = K
        ac(k) = sum(centered(1:N-k).*centered((k+1):N));
        ac(k) = ac(k)/(gamma0*(N-k));
    end
end