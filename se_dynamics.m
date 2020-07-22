% State estimation combined dynamics
% full state (x, x_hat, P)
%  x - true state
%  x_hat - state estimate
%  P - covariance

%% Example LQG System
dt = 0.2; N = 50; % trajectory discretization and length
A = 2*eye(2); B = dt*eye(2); C = eye(2); % process and measurement matrices
Q = diag([0.1 0.1]); % process noise covariance
R = diag([0.2 0.2]); % measurement noise covariance
K = dlqr(A,B,eye(2),eye(2)); % feedback law u = -Kx
x0 = [0;0]; P0 = eye(2); % initial state and covariance

u_nom = [1; 2]; % desired control

% iteratively generate nominal trajectory
x_nom = zeros(2,N); x_nom(:,1) = x0;
for i = 2:N
    x_nom(:,i) = A*x_nom(:,i-1) + B*u_nom;
end

%% Standard equations
x_est = x0; P = P0;
x = mvnrnd(x_est', P0, 1)';
X = nan(2,N);
Z = nan(2,N);
X_est = nan(2,N);
for i = 2:N
        
    % apply feedback control
    err = x_est - x_nom(:,i-1);
    u = u_nom - K*err;

    % dynamics
    w = mvnrnd([0 0], Q, 1)';
    x = A*x + B*u + w;
    X(:,i) = x;

    % noisy/biased measurement
    z = C*x + mvnrnd([0 0], R, 1)';
    Z(:,i) = z;

    % Kalman filter predict
    x_pred = A*x_est + B*u;
    P_pred = A*P*A' + Q;

    % Kalman filter update
    L = P_pred*C'/(C*P_pred*C' + R);
    x_est = x_pred + L*(z - C*x_pred);
    P = P_pred - L*C*P_pred;
    X_est(:,i) = x_est;
end

hold on
plot(X(1,:),X(2,:));
plot(Z(1,:),Z(2,:));
plot(X_est(1,:),X_est(2,:));

%% Combined dynamics equations
% for i = 2:N
%         
%     % apply feedback control
%     err = x_est - x_nom(:,i-1);
%     u = u_nom - K*err;
% 
%     % dynamics
%     w = mvnrnd([0 0], Q, 1)';
%     x = A*x + B*u + w;
%     X(:,i) = x;
% 
%     % noisy/biased measurement
%     z = C*x + mvnrnd([0 0], R, 1)';
% 
%     % Kalman filter predict
%     x_pred = A*x_est + B*u;
%     P_pred = A*P*A' + Q;
% 
%     % Kalman filter update
%     L = P_pred*C'/(C*P_pred*C' + R);
%     x_est = x_pred + L*(z - C*x_pred);
%     P = P_pred - L*C*P_pred;
% end

%% Continuous dynamics

C = eye(2); Q = eye(2); R = 2*eye(2);

xh0 = [0;0]; P0 = 5*eye(2);
%x0 = mvnrnd(xh0', P0, 1)';
x0 = [0;0];
s0 = [x0; xh0; vec(P0)];
u = [1;2];
tspan = [0 1];

[t_nom,x_nom] = ode45(@(t,x) dynamics(u),tspan,x0);

% [t,x] = ode45(@(t,x) dynamics(u,Q),tspan,x0);
% plot(x(:,1),x(:,2))

[t,s] = ode45(@(t,s) full_dynamics(s,u,C,Q,R),tspan,s0);

hold on
plot(s(:,1),s(:,2))
plot(s(:,3),s(:,4))

function x_dot = dynamics(u)
    x_dot = u;
end

function x_dot = noisy_dynamics(u,Q)
    w = mvnrnd([0;0],Q)';
    x_dot = u + w;
end

function u = control(t_nom,u_nom,t,)
    u = 
end

function s_dot = full_dynamics(s,u,C,Q,R)
    % extract state
    x = s(1:2);
    xh = s(3:4);
    P = reshape(s(5:end),2,2);
    
    % sample noise
    w = mvnrnd([0;0],Q)';
    v = mvnrnd([0;0],R)';
    
    % compute derivatives
    x_dot = u + w;
    L = P*C'*inv(R);
    xh_dot = u + L*(C*x + v - C*xh);
    P_dot = Q - L*R*L';
    s_dot = [x_dot; xh_dot; vec(P_dot)];
end