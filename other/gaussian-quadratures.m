% Computes Gauss quadrature points through the Legendre Polynomial roots. Also computes
% weights
% the idea of algorithm is pretty simple and it works out the legendre roots
% using iterative approach base on the number of points given. So this
% argument can be varied to obtain different weights and points in [-1,1].
% In principle, it can work for quite big N.

% Reference for algorithm:  Numerical Recipes in Fortran 77, Cornell press.

%set number of points
points = 3;
%initialize the points and the weights
x_p = zeros(points,1);                                           
weights = zeros(points,1);
nn = (points + 1)/2;
for i=1:nn
    initial_guess = cos(pi*(i-.25)/(points +.5));                        
    guess = initial_guess+1;
while abs(initial_guess - guess) > eps
    q0 = 1;
    q1 = 0;
    %compute legendre polynomials
    for j = 1:points
        q2 = q1;
        q1 = q0;
        q0 = ((2*j - 1)*initial_guess*q1-(j - 1)*q2)/j;              
    end
    p = points*(initial_guess*q0 - q1)/(initial_guess^2-1);
    guess = initial_guess;
    initial_guess = guess-q0/p;
end
    %collect points
    x_p(i) = -initial_guess;                                   
    x_p(points + 1 - i) = initial_guess;
    %collect weights
    weights(i) = 2/((1-initial_guess^2)*(p^2));                     
    weights(points + 1 - i) = weights(i);
end