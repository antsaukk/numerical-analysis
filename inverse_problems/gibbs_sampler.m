clear all;
close all;
% Gibbs sampler to reconstruct posterior probability distribution of the
% particle moving in the unit disk
%% Preallocations
f        = @(x,p) 1./abs(x-p);                                                   % charge potential function
post     = @(x,q,v,s,p) exp(-norm(v-q*f(x,p))^2/(2*s^2));                        % sampling posterior distribution

N        = 10000;                                                                % number of samples

q        = 1;                                                                    % charge of particle
x0       = 0.5 + 0.4i;                                                           % location of the particle
sigma    = 0.15;                                                                 % measurement variance
n        = 3;                                                                    % number of sensors
theta    = linspace(0, 2*pi, n+1);                                               % uniform intervals
p        = exp(1i*theta(1:end-1));                                               % sensor locations
sigma_n  = sigma*abs(max(f(x0,p)));                                              % noise in measurement
v        = q*f(x0,p) + sigma_n*randn(n, 1);                                      % measurement with noise

t        = 0:pi/50:2*pi;                                                         % points to generate boundary of unit disk

xgp      = -1:0.01:1;                                                            % x-grid points
ygp      = xgp;                                                                  % y-grid points
[X, Y]   = meshgrid(xgp, ygp);                                                   % grid
compl_pl = complex(X, Y);    
%% Gibbs Sampler
sample_histories = zeros(N, 2);                                                  % record sample histories
posterior        = zeros(size(compl_pl));                                        % value of posterior at every grid point

xk = [0, 0i];
for k = 1:N
    y = zeros(1, 2);
    for j = [1, 2]
        
    end
end
