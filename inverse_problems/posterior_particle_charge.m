clear all; 
close all;
%% Preallocations
q = 1; % charge of particle 
x = 0.5 + 0.4i; % true location of the particle
q0 = 1.1; % mean of the distribution of the charge
sigma_q = 0.2; % standard deviation of the distribution of the particle
%% (a), (b) 
% Posterior density when q is known is 
% \ksi_B1(x) exp^(1/(2 \sigma^2) \sum_{j=1}^{n} \norm(y - q f(x))^2)
% Posterior density when q is unknown 
% \ksi_B1(x) \sqrt(pi) exp^(B^2/(4A)) / \sqrt(A), A and B are defined in
% the code in detail
for n = [3, 9]
    theta = linspace(0, 2*pi, n);
    p = exp(i*theta);
    v = q./abs(x - p); % measurement
    for sigma = [0.05 0.15]
        sigma_n = sigma * abs(max(v));
        v_n = v + sigma_n * randn(size(v));
    end
end
