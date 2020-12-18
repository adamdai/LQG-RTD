% Simulate LQG trajectories 
% 
% Inputs:
%   sys    - struct containing system matrices A,B,and C, LQR gain K, and  
%            noise covariance matrices Q and R  
%   N      - trajectory length
%   N_traj - number of trajectories to generate
%   u_nom  - nominal input (constant)
%   x0     - initial state
%   P0     - initial covariance
%
% Outputs:
%   X - n x N x N_traj array containing N_traj n-dimensional trajectories 
%       of length N
%   

function X = simulate_LQG_trajectory(sys,N,N_traj,u_nom,x0,P0)
    
    % retrieve system matrices from sys struct
    A = sys.A; B = sys.B; C = sys.C;
    K = sys.K; Q = sys.Q; R = sys.R;
    
    n = size(A,1); % system dimension
    
    % iteratively generate nominal trajectory
    x_nom = zeros(n,N); x_nom(:,1) = x0;
    for i = 2:N
        x_nom(:,i) = A*x_nom(:,i-1) + B*u_nom;
    end
    
    X = nan(n,N,N_traj);

    % create N_traj trajectories
    for i_traj = 1:N_traj
    
        % initialize state and covariance
        x_est = x0; P = P0;
        x = mvnrnd(x_est', P0, 1)';

        % initialize trajectory
        X(:,1,i_traj) = x;

        for i = 2:N

            % apply feedback control
            err = x_est - x_nom(:,i-1);
            u = u_nom - K*err;

            % dynamics
            w = mvnrnd(zeros(n,1), Q, 1)';
            x = A*x + B*u + w;
            X(:,i,i_traj) = x;
            
            % noisy/biased measurement
            v = mvnrnd(zeros(n,1), R, 1)';
            z = C*x + v;

            % Kalman filter predict
            x_pred = A*x_est + B*u;
            P_pred = A*P*A' + Q;

            % Kalman filter update
            L = P_pred*C'/(C*P_pred*C' + R);
            x_est = x_pred + L*(z - C*x_pred);
            P = P_pred - L*C*P_pred;

        end
    end

end