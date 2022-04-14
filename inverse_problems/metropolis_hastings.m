clear all; 
close all;
% 
% Metropolis-Hasting Algorithm to reconstruct posterior
% probability density of particle moving inside the unit disk. 
%% Preallocations
F        = @(x) 1/(1-abs(x)^2)*x;                                                % unit disk mapping
F_inv    = @(x) (-1+sqrt(1+4*abs(x)^2))/(2*abs(x)^2)*x;                          % inverse of unit disk mapping
f        = @(x,p) 1./abs(x-p);                                                   % charge potential function
detp     = @(x,y) (1+abs(x)^2)*(1-abs(y)^2)^3/((1-abs(x)^2)^3*(1+abs(y)^2));     % first component of balance equation
expp     = @(x,y,q,v,s,p) exp((-norm(v-q*f(y,p))^2+norm(v-q*f(x,p))^2)/(2*s^2)); % second component of balance equation
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
compl_pl = complex(X, Y);                                                        % complex plane
%% Metropolis-Hastings Algorithm
iter = 1; 

% MH 
for gamma = [0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 1, 5]
    sample_histories = zeros(N, 2);                                  % record sample histories
    posterior        = zeros(size(compl_pl));                        % value of posterior at every grid point
    acceptance_rate  = 0;                                            % acceptance rate

    xk = 0 + 0i;                                                     % starting point of random walk
    for k = 1:N
        z   = F(xk);                                                 % Given x, compute F(x)
        w   = gamma * randn(1, 2);                                   % Draw W from N(0, gamma^2*I)
        wi  = complex(w(1), w(2));
        y   = F_inv(z + wi);                                         % Compute the next random step
        alp = expp(xk,y,q,v,sigma_n,p) * detp(xk,y);                 % Compute acceptance ratio
        tau = rand;                                                  % Draw from Unif(0,1)

        if tau < min(1,alp)                                          % Check the acceptance criteria
            xk              = y;
            acceptance_rate = acceptance_rate + 1;
        end
    
        sample_histories(k,1) = real(xk);                            % Real part of sample vector
        sample_histories(k,2) = imag(xk);                            % Imaginary part of sample
    end

    
    for i = 1:length(Y)                                              % compute posterior distribution
        for j = 1:length(X)
            KsiB           = double(abs(compl_pl(i,j)) < 1);         % characteristic function
            posterior(i,j) = KsiB*post(compl_pl(i,j),q,v,sigma_n,p); % posterior distribution 
        end
    end

    figure(iter)
    imagesc([-1,1], [-1,1], posterior);
    set(gca,'YDir','normal');
    axis image
    hold on
    plot(cos(t), sin(t), 'k');
    plot(p, '*');
    plot(x0, '+');
    title('Posterior');
    hold off
        
    iter = iter+1;

    figure(iter)
    hold on
    scatter(sample_histories(:,1), sample_histories(:,2),25,"red");
    plot(cos(t), sin(t), 'k');
    plot(p, '*', 'Color', 'b');
    plot(x0, '+', 'Color', 'k');
    cap = sprintf('Gamma=%.3f and acceptance rate=%.2f', gamma, acceptance_rate/N);
    title(cap, 'FontSize', 10);
    hold off

    iter = iter+1;

    figure(iter)
    hold on
    subplot(2,1,1);plot(sample_histories(:,1));
    title('Sample history x');
    subplot(2,1,2);plot(sample_histories(:,2));
    title('Sample history y');
    hold off

    iter = iter+1;

    ac_x = autocovariance(sample_histories(:,1));
    ac_y = autocovariance(sample_histories(:,2));
    L    = 100;

    figure(iter)
    hold on
    subplot(2,1,1);plot(0:L, ac_x(1:(L+1)));
    cap = sprintf('Autocovariance of x sample component with gamma=%.2f', gamma);
    title(cap, 'FontSize', 10);
    subplot(2,1,2);plot(0:L, ac_y(1:(L+1)));
    cap = sprintf('Autocovariance of y sample component with gamma=%.2f', gamma);
    title(cap, 'FontSize', 10);
    hold off

    iter = iter+1;
end

% We observe that acceptance rate is inversly proportional to sqrt(variance):
% gamma parameter. That is, when gamma is small, acceptance rate is high.
% On the other hand, when gamma is getting larger, acceptance rate
% decreases. One natural explanation is that small variance always
% induces small random step W, which never steps too far from the density.
% On the other hand, large variance results in highly dispersed steps and
% sampling points near the boundary is more likely, hence rejections are
% more likely, too. gamma from {0.05, 0.1, 0.3, 0.5, 1} demonstrates
% reasonable results and sampled observations resemble posterior density
% well. However, gamma from {0.001, 0.01} produces too dense reconstruction
% and gamma = 5 too sparse. These choices are probably not so good, since in
% the first case sampling does not really step away from initial point and
% resulting distribution is too dense, while in the second case too many
% points are sampled from the boundary which results in too low acceptance
% rate.
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