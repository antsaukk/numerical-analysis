%NEWTON - classical second order
f_x = @(x) (3*x-2)*(x*x+1);
dfx = @(x) 9*x*x - 4*x + 3;
x0 = 26; 
interval_length = 99; 
eps = 1e-8;
x(1) = x0; 
n = 2; 
iterations = interval_length + 1; 
while (n <= interval_length + 1)
  Fxn = f_x(x(n - 1));
  dFxn = dfx(x(n - 1));
  x(n) = x(n - 1) - Fxn/dFxn;
  if (abs(Fxn) <= eps)
    iterations = n; 
    break;
  end
  n = n + 1;
end
figure(1);
plot(0:iterations - 1,x(1:iterations),'bo-')
%hold on;
title('Convergence of the Newton method.')
legend('convergence of x_{n} approximations')
xlabel('iterations')
ylabel('x_{n}')
