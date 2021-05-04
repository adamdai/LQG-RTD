%% description
% Solving the uncertain normal equations A'*A*x = A'*b where A is an
% interval matrix and b is a vector (as standard)

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

b = [2; -1; 1];

%% sample solutions
% also sample inverses and eigenvalues
N = 10000; % number of A's to sample

X = zeros(2,N);
Ainv_lb = zeros(2,1);
Ainv_ub = zeros(2,1);
Aeig_set = zeros(2,1,N);

Aeig_lb = zeros(2,1);
Aeig_ub = zeros(2,1);
Aeig_set = zeros(2,1,N);

A_sample = randomSampling(A_z,N_A);

for i = 1:N
    A = A_sample{i};
    x = A\b;
    X(:,i) = x;
end

figure(); axis equal; grid on; hold on
scatter(reshape(X(1,:,:),[1,N_A*N_b]), reshape(X(2,:,:),[1,N_A*N_b]), '.')
xlabel('$x_1$','Interpreter','latex'); ylabel('$x_2$','Interpreter','latex');

