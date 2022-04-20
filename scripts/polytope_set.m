
b_c = [3; 2; 3; 2; 1];
b_g = rand(5,5);
b = zonotope(b_c, b_g);

A = [1 0 -1 0 1; 0 1 0 -1 1]';

center_poly = mptPolytope(A,b_c);

N_b = 1000;
b_sample = sampleBox(b,N_b);

hold on
plot(center_poly);
for i = 1:N_b
    b_i = b_sample(:,i);
    plot(mptPolytope(A,b_i));
end