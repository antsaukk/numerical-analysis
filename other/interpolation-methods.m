% (a)
%set the equispace grid = interpolation points
x = [-5:5/6:5]; 
%function to interpolate
f_x = 1./(1 + 25*x.^2); 
%plot function
figure(1);
plot(x, f_x, 'bo');
hold on;
%more dense grid to get better accuracy
x_0 = [-5:0.01:5];
%interpolating function on denser grid = more interpolation points
f_x1 = 1./(1 + 25*x_0.^2); 
%plot function on second grid
plot(x_0, f_x1, 'k');
%interpolation using polyfit and polyval
p = polyfit(x, f_x, 12); 
polynom = polyval(p, x_0); 
%plot the interpolation
plot(x_0, polynom, 'r');
title('Interpolating polynomial')

%(b) same task but with Chebyshev points
%Chebyshev points j = 0, ....,12 
j = [0:1:12]; 
x = 5*cos(pi*j/12); 
%Interpolating function to be computed on Chebyshev points
f_x = 1./(1 + 25*x.^2); 
figure(2)
plot(x, f_x, 'bo');
hold on;
%Same interpolating function as in (a)
f_x1 = 1./(1 + 25*x_0.^2); 
plot(x_0, f_x1, 'k');
%Interpolating polynomial
p = polyfit(x, f_x, 12); 
polynom = polyval(p, x_0); 
%plot the interpolation
plot(x_0, polynom, 'r');
title('Interpolating polynomial on Chebyshev points')

%(ñ) Barycentric coordinates to eliminate the ill-conditioning 
%Chebyshev points j = 0, ....,20 
j = [0:1:20]; 
x = 5*cos(pi*j/20); 
f_x = 1./(1 + 25*x.^2); 
figure(3);
plot(x, f_x, 'ko');
hold on;
%Berycentric weights
w = (-1).^j; 
w(1) = w(1)/2; 
w(21) = w(21)/2;
%Same interpolating function as in (a)
f_x1 = 1./(1 + 25*x_0.^2); 
plot(x_0, f_x1);
%Barycentric interpolation polynomial 
p = @(arg) sum(w.*f_x./(arg - x)) / sum(w./(arg - x));
%use arrayfun to evaluate function at x_0
bar_polynom = arrayfun(p, x_0); 
plot(x_0, bar_polynom, 'b');
title('Interpolating polynomial on Barycentric coordinates')
%Compute and plot the error
figure(4); 
error = abs(bar_polynom - f_x1); 
plot(x_0, error, 'r');
title('Interpolation error of barycentric coordinates')

%key observation there is that using Chebychev points gives
%significant improvement in interpolation estimate than using pure
%interpolation without them. However, even better result is obtained by
%increase of number of points and using barycentric coordinates.
%Barycentric coordinates could be useful when number of nodes is high and
%numerical software complains about ill-posed problem.
%Morale: Interpolation is is better to approach with Chebychev points since
%it helps to eliminate Runge's phenomenon, similar to the Gibbs phenomenon 
%that we can observe in Fourier series. For denser set of points one should
%use Barycentric. 


