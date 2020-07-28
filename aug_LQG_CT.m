% Augmented dynamics for continuous time LQG system
% state s:
%   [1:2]  x     - position
%   [3:4]  x_est - state estimates
%   [5:8]  P     - state estimate covariance
%   [9:10] k1,k2 - trajectory parameters
% 
% Inputs:
%   s - state
%
% Outputs:
%   ds - state derivative
%   

function ds = aug_LQG_CT(s,sys)

    % retrieve system matrices from sys struct
    A = sys.A; B = sys.B; C = sys.C;
    K = sys.K; Q = sys.Q; R = sys.R;
    
    % unpack state variable
    x = s(1:2);
    x_est = s(3:4);
    P = reshape(s(5:8),2,2);
    k = s(9:10);

    % get feedback control inputs
    u = k - K*(x_est - x_nom);

    % sample noise
    w = mvnrnd([0;0],Q)';
    v = mvnrnd([0;0],R)';

    % calculate the derivatives
    xd = A*x + B*u + w;
    L = P*C'/R; % Kalman gain
    x_estd = u + L*(C*x + v - C*x_est);
    Pd = Q - L*R*L';

    % return state derivative
    ds = [xd; x_estd; vec(Pd)] ;
    
end