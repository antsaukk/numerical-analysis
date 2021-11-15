clear all; 
close all;

h = 1/32;
t = 10; 
max = t / h; 

xy_v = zeros(max, 2); 
xy_v(1, :) = [1, 2]; 
xy_v1 = zeros(max, 2); 
xy_v1(1, :) = [2, 1]; 
time = zeros(max, 1); 

for i = 1:max
    xy_v(i+1, :) = RK4(xy_v(i, :), i*h, h);
    xy_v1(i+1, :) = RK4(xy_v1(i, :), i*h, h);
    time(i + 1) = time(i) + h; 
end

figure(1);
plot(time , xy_v(:, 1), 'o', time, xy_v(:, 2), 'x'); 
title('Populations development from IC (1, 2) versus time')
figure(2);
plot(xy_v(:, 1), xy_v(:, 2), 'o'); 
title('Phase portrain of the population in given predator-prey model')

figure(3);
plot(time , xy_v1(:, 1), 'o', time, xy_v1(:, 2), 'x'); 
title('Populations development from IC (2, 1) versus time')
figure(4);
plot(xy_v1(:, 1), xy_v1(:, 2), 'o');
title('Phase portrain of the population in given predator-prey model')



function next_iterate = RK4(x, t, h)
    %x
    K1 = f(x, t);
    %K1
    K2 = f(x + 0.5*h*K1', t + 0.5*h );
    %K2
    K3 = f(x + 0.5*h*K2', t + 0.5*h ); 
    %K3
    K4 = f(x + h*K3', t + h ); 
    %K4
    next_iterate = x + h/6 * (K1' + 2*K2' + 2*K3' + K4');
    %next_iterate
end
%first equation
function val=f(vals, t)
x1=vals(1); 
y1=vals(2); 
val=[2*x1*(1-x1/2 - y1/2);3*y1*(1 - 2*x1/3 - y1/3)];
end 

