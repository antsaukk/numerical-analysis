% Create 2D finite difference matrix
N = 50;
A = sparse(N^2,N^2);
h = 1/(N-1);
ijmap = @(i,j)((i-1)*N + j);
interior = []; % Non-boundary nodes

for i = 1:N
	for j = 1:N
        x(i, j) = (i-1)/(N-1);
        y(i,j) = (j-1)/(N-1);
        
        if ((i>1) & (i<N) & (j>1) & ( j<N))
            I1 = ijmap(i,j);
            interior = [interior I1];
            A(I1, ijmap(i-1,j)) =  -1/h^2;
            A(I1, ijmap(i+1,j)) =  -1/h^2;
            A(I1, ijmap(i,j-1)) =  -1/h^2;
            A(I1, ijmap(i,j+1)) =  -1/h^2;
            A(I1, I1)           =   4/h^2;
        end
    end
end

A = A(interior, interior);

L = chol(A, 'lower');
figure(1); subplot(2, 1, 1);
spy(L); title('Structure of Cholesky factor of FD matrix');
% Number of nonzeros
nnz(L);

order = amd(A);
L_amd = chol(A(order, order), 'lower');
subplot(2, 1, 2);
spy(L_amd); title('Structure of Cholesky factor of reordered FD matrix')
nnz(L_amd); % Returns 32911