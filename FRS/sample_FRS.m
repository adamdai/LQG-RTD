% Akshay's reachability formulation applied to a nominal trajectory 
% consisting of multiple linear segments

clc
clear
close all

%rng(252)

%% init variables

% % trajectory discretization and length
% t_f = 10; dt = 0.2; N = t_f/dt; 
% 
% % system
% A = eye(2); B = dt*eye(2); C = eye(2);
% K = dlqr(A,B,eye(2),eye(2)); % feedback law u = -Kx
% Q = diag([0.2 0.2]); % process noise covariance
% R = diag([0.1 0.1]); % measurement noise covariance

% trajectory discretization and length
t_f = 1; dt = 0.01; N = t_f/dt; 

% system
A = eye(2); B = dt*eye(2); C = eye(2);
K = dlqr(A,B,eye(2),eye(2)); % feedback law u = -Kx
Q = 1e-4*eye(2); % process noise covariance
R = 1e-4*eye(2); % measurement noise covariance

% form sys struct
sys.A = A; sys.B = B; sys.C = C;
sys.K = K; sys.Q = Q; sys.R = R;

x0 = [0;0]; % inital state
P0 = 1e-3*eye(2); % initial covariance

% trajectory parameter space
K1 = 2; K2 = 2; dk = 0.1;

% (uniform) sample nominal controls
N_samp = (2*K1/dk+1)*(2*K2/dk+1);
U_nom = cell(1,N_samp); i = 1;
for u1 = -K1:dk:K1
    for u2 = -K2:dk:K2
        U_nom{i} = repmat([u1; u2],1,N);
        i = i + 1;
    end
end


% generate nominal trajectories
X_nom = cell(1,N_samp);
for k = 1:N_samp
    u_nom = U_nom{k};
    x_nom = zeros(2,N); x_nom(:,1) = x0;
    for i = 2:N
        x_nom(:,i) = A*x_nom(:,i-1) + B*u_nom(:,i-1);
    end
    X_nom{k} = x_nom;
end

%%  run state estimator and generate reachable set along the way
% X = nan(2,N);
% Z = nan(2,N);
%     
% % initialize state and covariance
% x_est = x0; P = P0;
% x = mvnrnd(x_est', P0, 1)';
% 
% % initialize FRS
% pXrs = cell(1,M);
% 
% % for each segment
% for k = 1:M
%     
%     % compute FRS for this segment
%     if k == 1
%         X0 = probZonotope(x_est,cov2probGen(P),3);
%     else
%         X0 = pXrs{k-1}{N_seg};
%     end
%     %X0 = probZonotope(x_est,cov2probGen(P),3);
%     pXrs{k} = LQG_reach(X_nom{k},sys,X0);
%     
%     for j = 1:N_seg
%         
%         % full trajectory index
%         i = j + (k-1)*N_seg;
%         X(:,i) = x;
% 
%         % apply feedback control
%         err = x_est - x_nom(:,i);
%         u = u_nom(:,i) - K*err;
% 
%         % dynamics
%         w = mvnrnd([0 0], Q, 1)';
%         x = A*x + B*u + w;
% 
%         % noisy measurement
%         v = mvnrnd([0 0], R, 1)';
%         z = C*x + v;
%         Z(:,i) = z;
% 
%         % Kalman filter predict
%         x_pred = A*x_est + B*u;
%         P_pred = A*P*A' + Q;
% 
%         % Kalman filter update
%         L = P_pred*C'/(C*P_pred*C' + R);
%         x_est = x_pred + L*(z - C*x_pred);
%         P = P_pred - L*C*P_pred;
% 
%     end
% end


%% full trajectory reachability calculation

% union of all trajectory reachable sets at each time step
X_un = cell(1,N);

% initial reachable set X0
xZ0 = [0;0]; % (uncertain) mean
X0 = probZonotope(xZ0,cov2probGen(P0),3);

trajs = cell(1,N_samp);
for i = 1:N_samp
    trajs{i} = LQG_reach(X_nom{i},sys,X0);
end

%% plot 3-sigma confidence zonotopes
m = 3;
%Xrs = cell(1,N*N_samp);

figure(1)
axis equal;
hold on; grid on;

for i = 1:N
   for k = 1:N_samp
       umeanZ = mean(trajs{k}{i});
       covZ = cov2zonotope(sigma(trajs{k}{i}),m);
       Xrs = umeanZ + covZ;
%        if k == 1
%            X_un{i} = Xrs;
%        else
%            X_un{i} = convHull(X_un{i},Xrs);
%        end
       %X_un{i} = reduce(X_un{i},'girard',1);
       plt_FRS = plot(Xrs,[1,2],'Filled',true);
       plt_FRS.FaceAlpha = 0.1;
   end
%    plt_FRS = plot(X_un{i},[1,2],'Filled',true);
%    plt_FRS.FaceAlpha = 0.1;
end
