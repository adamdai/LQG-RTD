% Augmented dynamics for discrete time LQG system
% state s:
%   [1:2]   x     - position
%   [3:4]   x_est - state estimates
%   [5:8]   P     - state estimate covariance
%   [9:10]  k1,k2 - trajectory parameters
%   [11:12] x_nom - nominal position
% 
% Inputs:
%   ~ - dummy time
%   s - previous state
%   ~ - dummy control
%
% Outputs:
%   sn - next state
%   

function sn = aug_LQG_DT(tdummy,s,udummy,Tdummy)

    % system matrices
    dt = 0.2; A = eye(2); B = dt*eye(2); C = eye(2);
    K = dlqr(A,B,eye(2),eye(2)); Q = diag([0.1 0.1]); R = diag([0.1 0.1]);
    
    % unpack state variable
    x = s(1:2);
    x_est = s(3:4);
    P = reshape(s(5:8),2,2);
    k = s(9:10);
    x_nom = s(11:12);

    % get feedback control inputs
    u = k - K*(x_est - x_nom);

    % sample noise
    w = mvnrnd([0;0],Q)';
    v = mvnrnd([0;0],R)';
    
    % propagate forward nominal trajectory
    x_nom = A*x_nom + B*k;

    % dynamics and measurement model
    xn = A*x + B*u + w;
    y = C*xn + v;
    
    % Kalman predict
    x_estn = A*x_est + B*u;
    Pn = A*P*A' + Q;
    
    % Kalman update
    L = Pn*C'/(C*Pn*C' + R);
    x_estn = x_estn + L*(y - C*x_estn);
    Pn = Pn - L*C*Pn;
    
    % return combined next state
    sn = [xn; x_estn; vec(Pn); k; x_nom];
    
end