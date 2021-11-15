clear all; 
close all;

I = eye(3,3); 
T = [1 1 1];
eps = zeros(100,1);
k = zeros(100,1); 

k(1) = 2.0; 
for x = 2:length(k)
    k(x) = k(x-1) + 0.1; 
end


for i = 1:length(eps)
    eps(i) = 2^-k(i);
    
end

eps

err_one = zeros(100,1); %A - QR
err_two = zeros(100,1); %I - QQ

for l =1:length(eps) 
   A = [T; eps(l)*I]; 
   %A %check 

%modified GM 
Q = [];
R = [];
for i=1:size(A,2)
    q = A(:,i);
    for k=1:size(Q,2)
        R(k,i) = q'*Q(:,k);
        q = q - R(k,i)*Q(:,k);
    end
    R(i,i) = norm(q);
    Q(:,i) = q/R(i,i);
end

eA = abs(A-Q*R);    
err_one(l) = max(max(eA));
I_ = eye(3);
eI = abs(I_ - Q'*Q);
err_two(l) =  max(max(eI));    
end

figure(1);
loglog(eps, err_one); 
hold on; 
loglog(eps, err_two);
t = title('max||I-Q*QT|| and max||A-Q*R|| errors in modifed Gram-Schmidt');
t.FontSize = 9;
lg = legend('A-norm', 'I-norm');
lg.FontSize = 10;