%% description
% Solving the uncertain normal equations A'*A*x = A'*b where A is an
% interval matrix and b is a zonotope

%% setup
A_c = [1 6; 
       2 3;
       3 4];
%A_c = 5*rand(3,2);
%A_c = randi([-5,5],3,2);
A_g = [0.1 0.2; 
       0.2 0.5;
       0.1 0.4];
A_z = intervalMatrix(A_c,A_g);

b_c = [2; -1; 1];
b_g = 0.1*eye(3);
b_z = zonotope([b_c, b_g]);

%% sample solutions
N_b = 1000; % number of b's to sample
N_A = 1000; % number of A's to sample

X = zeros(2,N_A,N_b);

A_sample = randomSampling(A_z,N_A);
b_sample = sampleBox(b_z,N_b);

for i = 1:N_A
    for j = 1:N_b
        A = A_sample{i};
        b = b_sample(:,j);
        x = A\b;
        X(:,i,j) = x;
    end
end

figure(); axis equal; grid on; hold on
scatter(reshape(X(1,:,:),[1,N_A*N_b]), reshape(X(2,:,:),[1,N_A*N_b]), '.')
xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex');