% Evaluates nonlinear constraints for trajectory optimization
%  x - trajectory parameter choice
%  A_con - cell array of A matrix constraints
%  b_con - 

function [c, ceq, gc, gceq] = traj_opt_cons(x, c_k, g_k, A_con, b_con)
    epsilon = 1e-6;
    ceq = [];
    gceq = [];
    
    c = inf;
    lambdas = (x - c_k)./g_k; % given a parameter, get coefficients on k_slc_G generators
    for k = 1:length(A_con) % ============ since initial set is not sliceable w.r.t. inputs
        c_tmp = A_con{k}*lambdas - b_con{k}; % A*lambda - b <= 0 means inside unsafe set
        [c_tmp,i] = max(c_tmp); % max of this <= 0 means inside unsafe set
        c = min(c, c_tmp); % take smallest max. if it's <=0, then unsafe
        if c == c_tmp
            k_min = k;
            i_min = i;
        end
    end
    c = -c;
    
    gc = -A_con{k_min}(i_min,:)';
end