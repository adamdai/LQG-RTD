% simple distance-to-goal cost function for trajectory optimization 
% for linear system
%   sys   - struct containing system matrices
%   pos_i - position indices
%   N     - length of trajectory (number of timesteps)
%   k     - trajectory parameters
%   x_g   - desired goal
function [c, dc] = traj_opt_cost(sys,pos_i,N,k,x_g)

    % retrieve system matrices from sys struct
    A = sys.A; B = sys.B; 
    
    % calculate trajectory endpoint by propagating nominal dynamics
    x0 = zeros(size(A,1),1);
    D = 0;
    for i = 1:N
        D = D + A^(N-i)*B;
    end
    x = A^N * x0 + D * k; 
    
    % distance from trajectory endpoint to goal
    c = norm(x(pos_i) - x_g)^2;
%     dc = 2*[D(1,1)*(x(1) - x_g(1)) + D(2,1)*(x(2) - x_g(2)),
%             D(1,2)*(x(1) - x_g(1)) + D(2,2)*(x(2) - x_g(2))];
    dc = 2*D(pos_i,:)*(x(pos_i) - x_g);
end