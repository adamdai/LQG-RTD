clc
clear
close all

%rng(252)

%% init variables

% trajectory discretization and length
t_f = 10; dt = 0.2; N = t_f/dt; 

% segment time length and number of segments
t_seg = 2; M = t_f/t_seg; N_seg = t_seg/dt;

% system
A = eye(2); B = dt*eye(2); C = eye(2);
K = dlqr(A,B,eye(2),eye(2)); % feedback law u = -Kx
Q = diag([0.2 0.2]); % process noise covariance
R = diag([0.1 0.1]); % measurement noise covariance

% form sys struct
sys.A = A; sys.B = B; sys.C = C;
sys.K = K; sys.Q = Q; sys.R = R;

x0 = [0;0]; % inital state
P0 = eye(2); % initial covariance

% trajectory producing model 
% % constant control for time horizon
u_nom = [1; 1.5];
u_nom = repmat(u_nom,1,N);

% complex trajectory (composed linear segments)
%u_nom = [repmat([2;0],1,2*N_seg) repmat([1.5;1],1,2*N_seg) repmat([2;0],1,N_seg)];

% iteratively generate nominal trajectory
x_nom = zeros(2,N); x_nom(:,1) = x0;
for i = 2:N
    x_nom(:,i) = A*x_nom(:,i-1) + B*u_nom(:,i-1);
end

% divide nominal trajectory into segments
X_nom = cell(1,M);
for i = 1:M
    X_nom{i} = x_nom(:,N_seg*(i-1)+1:N_seg*i);
end

%%  run state estimator and generate reachable set along the way
X = nan(2,N);
Z = nan(2,N);
    
% initialize state and covariance
x_est = x0; P = P0;
x = mvnrnd(x_est', P0, 1)';

% initialize FRS
pXrs = cell(1,M);

% for each segment
for k = 1:M
    
    % compute FRS for this segment
    if k == 1
        X0 = probZonotope(x_est,cov2probGen(P),3);
    else
        X0 = pXrs{k-1}{N_seg};
    end
    %X0 = probZonotope(x_est,cov2probGen(P),3);
    pXrs{k} = LQG_reach(X_nom{k},sys,X0);
    
    for j = 1:N_seg
        
        % full trajectory index
        i = j + (k-1)*N_seg;
        X(:,i) = x;

        % apply feedback control
        err = x_est - x_nom(:,i);
        u = u_nom(:,i) - K*err;

        % dynamics
        w = mvnrnd([0 0], Q, 1)';
        x = A*x + B*u + w;

        % noisy measurement
        v = mvnrnd([0 0], R, 1)';
        z = C*x + v;
        Z(:,i) = z;

        % Kalman filter predict
        x_pred = A*x_est + B*u;
        P_pred = A*P*A' + Q;

        % Kalman filter update
        L = P_pred*C'/(C*P_pred*C' + R);
        x_est = x_pred + L*(z - C*x_pred);
        P = P_pred - L*C*P_pred;

    end
end


%% full trajectory reachability calculation

% initial reachable set X0
xZ0 = [0;0]; % (uncertain) mean
X0 = probZonotope(xZ0,cov2probGen(P0),3);

pXrs_full = LQG_reach(x_nom,sys,X0);

%% plot 3-sigma confidence zonotopes
m = 3;
Xrs_seg = cell(1,N);
Xrs_full = cell(1,N);

figure(1)
hold on; grid on;

% plot sampled trajectory
plot(X(1,:),X(2,:),'b');

% plot nominal trajectory
plot(x_nom(1,:),x_nom(2,:),'g');

% plot segment FRS
for k = 1:M % for each segment
    pXrs_seg = pXrs{k};
    color = rand(1,3);
    for j = 1:N_seg % plot zonotopes
        i = j + (k-1)*N_seg;
        umeanZ = zonotope(get(pXrs_seg{j},'Z'));
        covZ = cov2zonotope(get(pXrs_seg{j},'cov'),m);
        Xrs_seg{i} = umeanZ + covZ;
        plt_FRS = plotFilled(Xrs_seg{i},[1,2],color);
        plt_FRS.FaceAlpha = 0.2;
        plt_FRS.EdgeAlpha = 0.4;
    end
end

axis equal
xlabel('x-coordinate (m)');
ylabel('y-coordinate (m)');

figure(2)
hold on; grid on;

% plot sampled trajectory
plot(X(1,:),X(2,:),'b');

% plot nominal trajectory
plot(x_nom(1,:),x_nom(2,:),'g');

% plot full FRS
for i = 1:N
    umeanZ = zonotope(get(pXrs_full{i},'Z'));
    covZ = cov2zonotope(get(pXrs_full{i},'cov'),m);
    Xrs_full{i} = umeanZ + covZ;
    plt_FRS = plotFilled(Xrs_full{i},[1,2],[0.75,0.75,0.75]);
    plt_FRS.FaceAlpha = 0.2;
    plt_FRS.EdgeAlpha = 0.4;
end

axis equal
xlabel('x-coordinate (m)');
ylabel('y-coordinate (m)');

enclosed = 0;
for i = 1:N
    if in(Xrs_full{i},Xrs_seg{i})
        enclosed = enclosed + 1;
    end
end
fprintf("Percentage of Full FRS enclosed in segment-wise FRS: %f \n", enclosed / N);