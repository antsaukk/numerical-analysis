clear all;
close all;

%
[Q,R] = qr(rand(3)); 
A = Q*diag([10 0.5 0.1])*Q'

%(a) determine the eigenvectors and eigen values of A

[V,Diag] = eig(A)
V,Diag

%(b) compute A^ib as z0 = b, zi = Azi-1, 
%with b being a random
b = rand(3,1);
z_0 = b;
n=50;
alphas_1 = zeros(n,1);
alphas_2 = zeros(n,1);
alphas_3 = zeros(n,1);

for i = 1:n 
    z = A * z_0;
    alph = V\z;
    z_0 = z;
    alphas_1(i) = abs(alph(1));
    alphas_2(i) = abs(alph(2));
    alphas_3(i) = abs(alph(3));
    
end

figure(1);
clf;
%hold on;
semilogy(1:n, alphas_1)
hold on;
semilogy(1:n, alphas_2)
semilogy(1:n, alphas_3)
legend('alpha1', 'alpha2', 'alpha3')
    