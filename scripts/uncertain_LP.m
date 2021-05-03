%% parameters


%% Standard LP
% max_{x}     f^T x
% subject to  Ax <= b
%             x >= 0

f = [-1 -1/3];
A = [1 1
    1 1/4
    1 -1
    -1/4 -1
    -1 -1
    -1 1];

b = [2 1 2 1 -1 2];

x = linprog(f,A,b);

%% Uncertain LP

f = [-1 -1/3];
A = [1 1
    1 1/4
    1 -1
    -1/4 -1
    -1 -1
    -1 1];

b_c = [2 1 2 1 -1 2]';
b_g = 0.1*eye(6);
b_z = zonotope([b_c, b_g]);

N_b = 10000; % number of b's to sample

X = zeros(2,N_b);

options = optimoptions('linprog','Display','none');
b_sample = sampleBox(b_z,N_b);

for i = 1:N_b
    b_i = b_sample(:,i);
    x = linprog(f,A,b_i,[],[],[],[],options);
    X(:,i) = x;
end

scatter(X(1,:), X(2,:), '.')