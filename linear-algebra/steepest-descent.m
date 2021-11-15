clear all;
close all;

%steepest descent method to approximate Ax=b
A = [2 1; 1 5]; 
eps = 10^-12;
b = [1 1]'; 
solution = zeros(2,1);
r = b - A*solution;
p = r;
errors = [];
iter = 0;


for i=1:50
    solution
    errors(end + 1) = norm(r);
    iter = i;
    alpha = (p'*r)/(p'*A*p);
    solution = solution + alpha*p; 
    r = b - A*solution;
    p = r;
    if (norm(r) < eps)
        fprintf('converged to %d iterations\n', i)
        break;
    end
end

fprintf('steepest gradient solution');
solution
fprintf('backslash  solution');
A\b
semilogy(1:iter, errors, 'o--');
t = title('Errors in iteration of steepest descent');
t.FontSize = 9;