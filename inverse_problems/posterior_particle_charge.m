clear all; 
close all;
%% Preallocations
q = 1; % charge of particle 
x = 0.5 + 0.4i; % true location of the particle
q0 = 1.1; % mean of the distribution of the charge
sigma_q = 0.2; % standard deviation of the distribution of the particle

% grid of points to evaluate the posterior
xgp = -1:0.01:1;
ygp = xgp;
[X, Y] = meshgrid(xgp, ygp);
complex_plane = complex(X, Y);

% unit disk visualizaion
t = 0:pi/50:2*pi;

% store values of both posteriors
Posterior_a = zeros(size(complex_plane));
Posterior_b = zeros(size(complex_plane));
%% (a), (b) 
% Posterior density when q is known is 
% \ksi_B1(x) exp^(1/(2 \sigma^2) \sum_{j=1}^{n} \norm(y - q f(x))^2)
% Posterior density when q is unknown 
% \ksi_B1(x) \sqrt(pi) exp^(B^2/(4A)) / \sqrt(A), A and B are defined in
% the code in detail
iter = 1;
for n = [3, 9]
    % uniform intervals
    theta = linspace(0, 2*pi, n+1);

    % sensor locations
    p = exp(1i*theta(1:end-1));
    
    % measurement without noise
    v = q./abs(x - p); 
    for sigma = [0.05 0.15]
        % noise
        sigma_n = sigma * abs(max(v)); 

        % noised measurements
        v_n = v + sigma_n * randn(1,n);

        % evaluation of posteriors
        for i = 1:length(Y)
            for j = 1:length(X)
                % characteristic function
                KsiB = double(abs(complex_plane(i,j)) < 1);

                % f helper
                 f = 1./abs(complex_plane(i,j) - p);

                % posterior distribution for (a)
                Posterior_a(i,j) = KsiB*exp(-1/(2*sigma_n^2)*norm(v_n-q*f)^2);

                % posterior distribution for (b)
                A  = 1/(2*(sigma_n*sigma_q)^2)*(sigma_q^2*norm(f)^2+sigma_n^2);
                B  = 1/(sigma_n*sigma_q)^2*(sigma_q^2*v_n*f' + sigma_n^2*q0);
                Posterior_b(i,j) = KsiB*sqrt(pi)*exp(B^2/(4*A))/sqrt(A);
            end
        end

        % posterior a visualization
        figure(iter)
        imagesc([-1,1], [-1,1], Posterior_a);
        axis image
        hold on
        plot(cos(t), sin(t), 'k');
        plot(p, '*');
        plot(x, '+');
        caption1 = sprintf('Posterior (a) with n=%d and std=%.3f', n, sigma_n);
        title(caption1, 'FontSize', 10);
        hold off

        iter=iter+1;

        % posterior b visualization
        figure(iter)
        imagesc([-1,1], [-1,1], Posterior_b);
        axis image
        hold on
        plot(cos(t), sin(t), 'k');
        plot(p, '*');
        plot(x, '+');
        caption2 = sprintf('Posterior (b) with n=%d, std=%.3f and E[q]=%.3f, std(q)=%.3f', n, sigma_n, q0, sigma_q);
        title(caption2, 'FontSize', 10);
        hold off

        iter=iter+1;
    end
end

% Observe that when variance is small and there are 9 positioned sensors,
% they are able to identify the location of the charge fairly accurately
% and posterior probability function is densely concentrated near the actual
% location. When the variance grows, 9 sensors are still conduct pretty
% accurate measurements, however the density function is more dispearsed. 
% Relatively same behaviour is observed for 3 sensors, however, the accuracy
% is now clearly worse, which is reasonable, since there are less
% information available.
