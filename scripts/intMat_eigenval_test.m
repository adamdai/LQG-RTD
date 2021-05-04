% compare overapproximated interval matrix inverse and true inverse set

C = [2 3; 4 5];
G = [1 0.5; 0.2 2];
A = intervalMatrix(C,G);

A = intMat_transpose(A) * A;
   
% sample for intMat
N = 10000;
A_samp = randomSampling(A,N);
Aeig_lb = zeros(2,1);
Aeig_ub = zeros(2,1);
Aeig_set = zeros(2,1,N);
for i = 1:N
    Ai = A_samp{i};
    Ai_eig = eig(Ai);
    % update lb
    Aeig_lb = min(Aeig_lb, Ai_eig);
    % update ub
    Aeig_ub = max(Aeig_ub, Ai_eig);
    % store points
    Aeig_set(:,:,i) = Ai_eig;
end

% plot overapproximate intervals for each eigenvalue along with true sample
% points

for i = 1:2
    subplot(2,1,i); hold on
    plot([Aeig_lb(i);Aeig_ub(i)],[0;0],'-o','LineWidth',2);
    xlim([Aeig_lb(i)-1, Aeig_ub(i)+1]); ylim([-1 1]);
    scatter(Aeig_set(i,:),zeros(N,1),'*');
    title(strcat('\lambda_',num2str(i)));
end







