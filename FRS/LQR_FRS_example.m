clc
clear
close all

rng(252)

%% init variables

% trajectory discretization and length
t_f = 10; dt = 0.2; N = t_f/dt; 

% segment time length and number of segments
t_seg = 2; M = tf/t_seg; N_seg = t_seg/dt;

% system
A = eye(2); B = dt*eye(2); C = eye(2);
K = dlqr(A,B,eye(2),eye(2)); % feedback law u = -Kx
Q = diag([0.1 0.4]); % process noise covariance
R = diag([0.2 0.1]); % measurement noise covariance

% form sys struct
sys.A = A; sys.B = B; sys.C = C;
sys.K = K; sys.Q = Q; sys.R = R;

x0 = [0;0]; % inital state
P0 = eye(2); % initial covariance

% define trajectory parameter space
K1 = 0.2; % -0.2 to +0.2
K2 = 0.2; % -0.2 to +0.2

% compute nominal trajectory FRS
x_nom = linear_FRS(K1,K2,t_f,dt);

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
        plotFilled(Xrs_seg{i},[1,2],color);
    end
end

% plot full FRS
for i = 1:N
    umeanZ = zonotope(get(pXrs_full{i},'Z'));
    covZ = cov2zonotope(get(pXrs_full{i},'cov'),m);
    Xrs_full{i} = umeanZ + covZ;
    plotFilled(Xrs_full{i},[1,2],[0.75,0.75,0.75]);
end

enclosed = 1;
for i = 1:N
    if ~in(Xrs_full{1},Xrs_seg{1})
        enclosed = 0;
    end
end
fprintf("Full FRS enclosed in segment-wise FRS: %d \n", enclosed);

% axis equal
% xlabel('x-coordinate (m)');
% ylabel('y-coordinate (m)');
% 
% figure(2)
% hold on; grid on;
% 
% % plot sampled trajectory
% plot(X(1,:),X(2,:),'b');
% 
% % plot nominal trajectory
% plot(x_nom(1,:),x_nom(2,:),'g');

axis equal
xlabel('x-coordinate (m)');
ylabel('y-coordinate (m)');