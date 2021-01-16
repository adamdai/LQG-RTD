% simple distance-to-goal cost function for trajectory optimization 
% for linear system
%   sys   - struct containing system matrices
%   pos_i - position indices
%   N     - length of trajectory (number of timesteps)
%   k     - trajectory parameters
%   x_g   - desired goal
%   M     - trajectory parameter to position map from linear planning model
function [c, dc] = traj_opt_cost(sys,pos_i,N,k,ic,x_g,M)

    % retrieve system matrices from sys struct
    A = sys.A; B = sys.B; 
    
    k = [ic(1); ic(3); k(1); ic(2); ic(4); k(2)];
    
    % calculate trajectory endpoint by propagating nominal dynamics
    x0 = zeros(size(A,1),1);
    D = 0;
    for i = 1:N
        D = D + A^(N-i) * B * M{i};
    end
    x = A^N * x0 + D * k; 
    
    % distance from trajectory endpoint to goal
    c = norm(x(pos_i) - x_g)^2;
    dc = 2*[D(1,3)*(x(1) - x_g(1)) + D(2,3)*(x(2) - x_g(2)),
            D(1,6)*(x(1) - x_g(1)) + D(2,6)*(x(2) - x_g(2))];
    %dc = 2*D(pos_i,:)*(x(pos_i) - x_g);
end