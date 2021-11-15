clear all; 
close all; 
M = transpose([1 1 0; 1 0 1]); 
b = transpose([1 0 0]); 
M
b
minv = M\b; 
fprintf("comparison to the results computed by hand"); 
minv
%(a) Find an orthonormal basis for R(A) by using Gram-Schmidt process
%(b) Find the QR-decomposition of A 
fprintf("TASK4"); 
A = [1 2 2; 
     1 1 2; 
     2 1 1; 
     2 2 1]; 
A
my_gsmith(A); 
function[Q,R] = my_gsmith(A) 
Q = [];
for i=1:size(A,2)
    q = A(:,i);
    %q is the column
    for k=1:size(Q,2)
        R(k,i) = q'*Q(:,k);
        q = q - R(k,i)*Q(:,k);
    end
   R(i,i) = norm(q);
   u = q/R(i,i);
   fprintf("the i'th=%d column vector is", i);  
   u
   Q(:,i) = q/R(i,i); 
end
fprintf("the [Q,R] decomposition is");
Q
R 
%test function
fprintf("conduct the normalization tests\n");
for y=1:size(Q, 2)
    q1 = Q(:,y); 
    s = 0;
    for x = 1:size(q1)
        s = s + q1(x)^2; 
    end
    fprintf("square root of the y'th=%d column in q equals %d\n", y, sqrt(s));
end
%S = sqrt(0.3162^2 +  0.3162^2 + 0.6325^2 + 0.6325^2);
%S
end
 



    