clear all; 
close all; 

%(d) solving with FE - modified
h = 1/1000; 
u = [-sqrt(2), 17]; 
t = 4/h; 
FEm = zeros(t, 2); 
time = zeros(t, 1); 
%for n = 1:3
    for i = 1:t
        y = (1 - h)^i .* u; 
        FEm(i,:) = y;
        time(i) = i*h; 
    end
    plot(time, FEm); 
    hold on;
    %plot(1:t, exact(:,1)); 
    %h = h*(1/10);
    %t = 4/h
%end

title('FE - modifed');
%legend('h = 1/100', 'h = 1/100', 'h = 1/1000', 'h = 1/1000', 'h = 1/10000', 'h = 1/10000');
legend('y1', 'y2');

%(d) solving with FE 
A = [16 sqrt(2); 0 -1];

F = @(y, h) (eye(2) + h .* A) * y;

h = 1/1000; 
y = [-sqrt(2), 17]'; 
t = 4/h; 
FE = zeros(t, 2);
time1 = zeros(t, 1); 
%for n = 1:3
    for i = 1:t
        y = F(y, h); 
        %y
        FE(i,:) = y;
        time1(i) = i*h; 
        %break
    end
    figure(2)
    plot(time1, FE);
    title('FE');
    legend('y1', 'y2');
    %plot(1:t, exact(:,1)); 
    %h = h*(1/10);
    %t = 4/h
%end

%(d) solving with RK4

f = @(x, t) [16 sqrt(2); 0 -1]*x';
h = 1/2000; 
t1 = 4 / h;
rk4 = zeros(t1, 2); 
y = [-sqrt(2), 17]; 
time2 = zeros(t, 1); 
for i = 1:t1
    y = RK4(f, y, i, h);
    rk4(i,:) = y;
    time2(i) = i*h; 
    %rk4(i,:)
end
figure(3)
plot(time2, rk4);
title('RK4');
legend('y1', 'y2');


function next_iterate = RK4(f, x, t, h)
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

