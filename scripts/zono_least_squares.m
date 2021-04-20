%% Standard LS

A = [2 1; 
     0 3;
     1 4];
b = [2; -1; 1];

x = inv(A'*A) * A'*b;


%% Zonotope LS (uncertain b)

A_z = [2 1; 
       0 3;
       1 4];
b_c = [2; -1; 1];
b_g = [  0   0;
       0.1   0;
         0 0.1];
b_z = zonotope([b_c, b_g]);

x_z = inv(A_z'*A_z) * A_z'*b_z;


