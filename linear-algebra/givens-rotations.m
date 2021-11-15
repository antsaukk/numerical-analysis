clear all; 
close all;

%steepest descent method
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

%Givens_rotations 
[Q, R] = my_givens_qr(A);    

eA = abs(A-Q*R); %is this how it should be?  
err_one(l) = max(max(eA));
I_ = eye(4); %or eye(4)
eI = abs(I_ - Q'*Q); %and this one??? 
err_two(l) =  max(max(eI));    
end

figure(1);
loglog(eps, err_one); 
hold on; 
loglog(eps, err_two);
t = title('max||I-Q*QT|| and max||A-Q*R|| errors in Givens Rotations');
t.FontSize = 9;
lg = legend('A-norm', 'I-norm');
lg.FontSize = 10;

function [Q,A] = my_givens_qr(A) 
   Q = eye(max(size(A)));
        for i=1:size(A,2)
            for j=(i+1):size(A,1)
                % Construct G
                x  = [A(i,i) ; A(j,i)];
                xN = [-x(2)  ; x(1)];
                G = [ x'/norm(x) ; xN'/norm(xN)];
                % Operate with G
                Q([i j],:) = G*Q([i j],:);
                A([i j],:) = G*A([i j],:);        
            end
        end
   Q = Q';
end
