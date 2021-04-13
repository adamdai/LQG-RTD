% returns the inverse of a matrix zonotope

% function Ainv = matZ_inverse(A)
%     
% end

C = [2 3; 4 5];
G = [1 0.5; 0.2 2];
A = intervalMatrix(C,G);

% sample for intMat
N = 1000;
A_samp = randomSampling(A,N);
Ainv_samp = {};
Ainv_lb = zeros(2,2);
Ainv_ub = zeros(2,2);
for i = 1:N
    Ai = A_samp{i};
    Ai_inv = inv(Ai);
    % update lb
    Ainv_lb = min(Ainv_lb, Ai_inv);
    % update ub
    Ainv_ub = max(Ainv_ub, Ai_inv);
end



