% Discrete time Turtlebot
% - discrete LQR
% - Kalman filter
% state: s = (x,y,h)
% input: u = (v,w)
% dynamics: (x_k,y_k,h_k) = (x_{k-1} + v*cos(h_{k-1}),y_{k-1} +
% v*sin(h_{k-1}),w) + N(0,Q)
% measurement model: (x_m, y_m) = (x,y) + N(0,R)


%% parameters

n = 3;                  % state dimension
m = 4;                  % measurement dimension
dt = 0.1;               % time discretization
Q = 0.1*dt*eye(n);      % process covariance
R = 0.1*eye(m);         % measurement covariance
v = 1;                  % robot scalar velocity

T = 0:dt:10;            % time interval
x0 = zeros(n,1);        % initial state
mu0 = zeros(n,1);       % initial state mean
Sigma0 = 0.01*eye(n);   % initial state covariance

M = [0 10 10  0;        % measurement beacon locations
     0  0 10 10];

%% Filter
x = x0;
xlist = [];
ylist = [];
mu = mu0;
Sigma = Sigma0;
mulist = [];

for t = T
    % Dynamics update
    u = [v; sin(t)];
    x = f(x,u,dt) + mvnrnd(zeros(n,1),Q)';
    y = g(x,M) + mvnrnd(zeros(m,1),R)';
    xlist = [xlist x];
    ylist = [ylist y];
    
    % EKF predict
    mu = f(mu,u,dt);
    % Jacobians
    At = A(mu,v,dt);
    Ct = C(mu,M);
    Sigma = At*Sigma*At' + Q;
    % EKF update
    mu = mu + Sigma*Ct'*inv(Ct*Sigma*Ct' + R)*(y - g(mu,M));
    Sigma = Sigma - Sigma*Ct'*inv(Ct*Sigma*Ct' + R)*Ct*Sigma;
    mulist = [mulist mu];

end

figure();
plot(T,xlist(1,:),T,mulist(1,:));
title('p^x');
legend('true','estimated');
figure();
plot(T,xlist(2,:),T,mulist(2,:));
title('p^y');
legend('true','estimated');
figure();
plot(T,xlist(3,:),T,mulist(3,:));
title('\theta');
legend('true','estimated');


%% functions

% dynamics
% x = [px; py; theta]
% u = [v; phi]
function x = f(x,u,dt)
    x = x + dt*[u(1)*cos(x(3)); u(1)*sin(x(3)); u(2)];
end

% measurement
function y = g(x,M)
    y  = vecnorm(M - x(1:2), 2, 1)';
end

% dynamics jacobian
function At = A(x,v,dt)
    At = [1 0 -dt*v*sin(x(3)); 
          0 1 dt*v*cos(x(3)); 
          0 0 1];
end

% measurement jacobian
function Ct = C(x,M)
    m = size(M,2); n = size(x,1);
    Ct = zeros(m,n);
    for i = 1:m
        d = norm(M(:,i) - x(1:2));
        Ct(i,:) = [-(M(1,i)-x(1))/d, -(M(2,i)-x(2))/d, 0];
    end
end