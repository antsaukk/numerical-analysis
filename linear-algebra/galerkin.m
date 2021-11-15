clear all;
close all;

%5.4
%(a) determine eigenvalues and eigenvectors of A 
n=20; 
L = 1000;  
%[Q,R] = qr(rand(n));
%A = Q*diag(linspace(1,L,n)*Q')
%[V,D] = eig(A)
b = rand(n,1)
errors = zeros(10,1);



i = 1
for k = 5:15
    [Qs,R] = qr(rand(n));
    A = Qs*diag(linspace(1,L,n))*Qs'
    Q = my_arnoldi(A,b,k)
    Xn = (Q'*A*Q)\(Q'*b);
    X = Q*Xn;
    error = X - A\b; 
    error_A = sqrt(error'*A*error)
    errors(i) = error_A 
    i = i + 1
end
errors
figure(1);
clf;
semilogy(5:15, errors, 'o--')


%funcion to compute Basis for krylov subspace
function [Q,R] = my_arnoldi(A,b,n)
    Q = [];
    q = b;
    for i=1:n
        for k=1:size(Q,2)
            R(k,i) = q'*Q(:,k);
            q = q - R(k,i)*Q(:,k);
        end
        R(i,i) = norm(q);
        Q(:,i) = q/R(i,i);
        q = A*Q(:,i);
    end
end