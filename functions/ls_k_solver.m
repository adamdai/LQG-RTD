% computes constrained least-squares solution to optimization of k using
% distance-to-goal cost function with no obstacle constraints
%   sys   - struct containing system matrices
%   pos_i - position indices
%   N     - length of trajectory (number of timesteps)
%   k     - trajectory parameters
%   x_g   - desired goal
%   M     - trajectory parameter to position map from linear planning model
function k = ls_k_solver(sys, pos_i, N, ic, x_g, M, lb, ub)

    % retrieve system matrices from sys struct
    A = sys.A; B = sys.B; 
    
    % form full D matrix
    D = 0;
    for i = 1:N
        D = D + A^(N-i) * B * M{i};
    end
    
    % split into D_peak and D_init
    D_peak = D(pos_i,[3 6]);
    D_init = D(pos_i,[1 2 4 5]);
    
    % least squares solution
    A_ls = D_peak;
    b_ls = x_g - D_init * ic; 
    k = lsqlin(A_ls, b_ls, [], [], [], [], lb, ub);
end