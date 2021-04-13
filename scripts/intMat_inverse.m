% returns the inverse of a matrix zonotope
% (cop out with sampling for now)

% C = [2 3; 4 5];
% G = [1 0.5; 0.2 2];
% A = intervalMatrix(C,G);

function Ainv = intMat_inverse(A)
   
% check for nonzero volume
if volume(A) == 0
    Ainv = intervalMatrix(inv(center(A)));
else
    % sample for intMat
    N = 1000;
    A_samp = randomSampling(A,N);
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

    Ainv_c = (Ainv_lb + Ainv_ub) / 2;
    Ainv_r = (Ainv_ub - Ainv_lb) / 2;

    Ainv = intervalMatrix(Ainv_c, Ainv_r);
end

end







