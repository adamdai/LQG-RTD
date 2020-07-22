clc
clear
% close all

%% init variables

dx = 1; dt = 0.2; N = 50; % trajectory discretization and length
A = eye(2); B = dt*eye(2); C = eye(2);
% tracking controller
K = dlqr(A,B,eye(2),eye(2)); % feedback law u = -Kx

x0 = [0;0]; % inital state
P0 = eye(2); % initial covariance
Q = 0.1*eye(2); % process noise covariance
R = 0.1*eye(2); % measurement noise covariance

% trajectory producing model 
% given x0, u_nom = (u1,u2), generate x_nom
u_nom = [1; 2]; % constant control for time horizon

% iteratively generate nominal trajectory
x_nom = zeros(2,N); x_nom(:,1) = x0;
for i = 2:N
    x_nom(:,i) = A*x_nom(:,i-1) + B*u_nom;
end

N_traj = 100; % number of sample trajectories

%% generate sample trajectories
X_tr = nan(2,N,N_traj); % tracking
X_se = nan(2,N,N_traj); % state estimating
X_est = nan(2,N,N_traj); % estimate

Z = nan(2,N,N_traj); % measurements

for i_traj = 1:N_traj
    
    % initialize state and covariance
    x_est = x0; P = P0;
    %x = mvnrnd(x_est', P0, 1)';
    x = x0;
    x_tr = x; x_se = x;
    
    % initialize trajectory
    X_tr(:,1,i_traj) = x_tr;
    X_se(:,1,i_traj) = x_se;
    X_est(:,1,i_traj) = x_est;
    
    for i = 2:N
           
        % sample process noise
        w = mvnrnd([0 0], Q, 1)';
        
        % tracking trajectory
        % apply feedback control
        err = x_tr - x_nom(:,i-1);
        %u = u_nom - K*err;
        u = u_nom;
        
        % dynamics
        x_tr = A*x_tr + B*u + w;
        X_tr(:,i,i_traj) = x_tr;
        
        % state estimating trajectory
        % apply feedback control
        err = x_est - x_nom(:,i-1);
        u = u_nom - K*err;
        
        % dynamics
        x_se = A*x_se + B*u + w;
        X_se(:,i,i_traj) = x_se;
        
        % noisy/biased measurement
        z = C*x_se + mvnrnd([0 0], R, 1)';
        Z(:,i,i_traj) = z;

        % Kalman filter predict
        x_pred = A*x_est + B*u;
        P_pred = A*P*A' + Q;

        % Kalman filter update
        L = P_pred*C'/(C*P_pred*C' + R);
        x_est = x_pred + L*(z - C*x_pred);
        P = P_pred - L*C*P_pred;
        
        X_est(:,i,i_traj) = x_est;
        
    end
    
end

%% plot trajectories

figure()
hold on; grid on;

% plot sampled trajectories
for i_traj = 1:N_traj
    plot(X_tr(1,:,i_traj),X_tr(2,:,i_traj),'b');
    plot(X_se(1,:,i_traj),X_se(2,:,i_traj),'r');
end

% plot nominal trajectory
plot(x_nom(1,:),x_nom(2,:),'g');

% % plot state estimates
% plot(X_est(1,:),X_est(2,:),'c');

fprintf("Tracking error (x - x_nom): %f \n", sum(vecnorm(X_tr - x_nom,2,1),'all')/N_traj);
fprintf("State estimation error (x_hat - x): %f \n", sum(vecnorm(X_est - X_se,2,1),'all')/N_traj);
fprintf("State estimating tracking error (x - x_nom): %f \n", sum(vecnorm(X_se - x_nom,2,1),'all')/N_traj);


%%
w_avg = 0;
Q = [2 1; 1 3];
N = 1e5;
for i = 1:N
    w = mvnrnd([0 0], Q, 1)';
    w_avg = w_avg + norm(w)^2;
end
disp(w_avg/N)

%%
S = zeros(2,50,100);
Q = eye(2);
for n = 1:1000
    s = zeros(2,1);
    for i = 1:50
        w = mvnrnd([0 0], Q, 1)';
        s = s + w;
        S(:,i,n) = s;
    end
end
fprintf("error: %f \n", sum(vecnorm(S,2,1),'all')/1000);