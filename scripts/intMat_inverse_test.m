% compare overapproximated interval matrix inverse and true inverse set

C = [2 3; 4 5];
G = [1 0.5; 0.2 2];
A = intervalMatrix(C,G);

A = intMat_transpose(A) * A;
   
% sample for intMat
N = 10000;
A_samp = randomSampling(A,N);
Ainv_lb = zeros(2,2);
Ainv_ub = zeros(2,2);
Ainv_set = zeros(2,2,N);
for i = 1:N
    Ai = A_samp{i};
    Ai_inv = inv(Ai);
    % update lb
    Ainv_lb = min(Ainv_lb, Ai_inv);
    % update ub
    Ainv_ub = max(Ainv_ub, Ai_inv);
    % store points
    Ainv_set(:,:,i) = Ai_inv;
end

% plot overapproximate intervals for each entry along with true sample
% points

% Ainv_(1,1)
for i = 1:2
    for j = 1:2
        subplot(4,1,2*i+j-2); hold on
        plot([Ainv_lb(i,j);Ainv_ub(i,j)],[0;0],'-o','LineWidth',2);
        xlim([Ainv_lb(i,j)-1, Ainv_ub(i,j)+1]); ylim([-1 1]);
        scatter(Ainv_set(i,j,:),zeros(N,1),'*');
        title(sprintf('[%d,%d]',i,j));
    end
end
% % Ainv_(1,2)
% subplot(4,1,2); hold on
% plot([Ainv_lb(1,1);Ainv_ub(1,1)],[0;0],'-o','LineWidth',2);
% xlim([Ainv_lb(1,1)-1, Ainv_ub(1,1)+1]); ylim([-1 1]);
% scatter(Ainv_set(1,1,:),zeros(N,1),'*');
% 
% % Ainv_(2,1)
% subplot(4,1,3); hold on
% plot([Ainv_lb(1,1);Ainv_ub(1,1)],[0;0],'-o','LineWidth',2);
% xlim([Ainv_lb(1,1)-1, Ainv_ub(1,1)+1]); ylim([-1 1]);
% scatter(Ainv_set(1,1,:),zeros(N,1),'*');
% 
% % Ainv_(2,2)
% subplot(4,1,4); hold on
% plot([Ainv_lb(1,1);Ainv_ub(1,1)],[0;0],'-o','LineWidth',2);
% xlim([Ainv_lb(1,1)-1, Ainv_ub(1,1)+1]); ylim([-1 1]);
% scatter(Ainv_set(1,1,:),zeros(N,1),'*');









