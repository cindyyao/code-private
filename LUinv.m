function pA = LUinv(A)

% perform LU factorization using matlab function
 [L, U]=lu(A);

 % solving Ax=[1; 0; 0], Ax=[0; 1; 0], and Ax=[0; 0; 1]
 d=L\[1;0;0];
 col1 = U\d;
 d=L\[0;1;0];
 col2 = U\d;
 d=L\[0;0;1];
 col3 = U\d;
 pA = [col1 col2 col3];
end