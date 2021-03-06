% example for computing FRS for LQG agent over a space trajectory
% parameters and intersecting it with a zonotopic obstacle to determine
% unsafe trajectory parameters

clc
clear
close all

%% 2D system
% n = 2; % dimension
% 
% % trajectory discretization and length
% t_f = 10; dt = 0.2; N = t_f/dt; 
% 
% % system
% A = eye(2); B = dt*eye(2); C = eye(2);
% K = dlqr(A,B,eye(2),eye(2)); % feedback law u = -Kx
% Q = diag([0.1 0.1]); % process noise covariance
% R = diag([0.1 0.1]); % measurement noise covariance
% 
% % form sys struct
% sys.A = A; sys.B = B; sys.C = C;
% sys.K = K; sys.Q = Q; sys.R = R;
% 
% x0 = [0;0]; % inital state
% P0 = 0.3*eye(2); % initial covariance
% 
% % define trajectory parameter space
% K1 = 5; 
% K2 = 4; 
% U_nom = zonotope([0 K1 0; 0 0 K2]);

%% double integrator system
n = 4; % state dimension

% trajectory discretization and length
t_f = 5; dt = 0.05; N = t_f/dt; 

% system
A = eye(n);
A(1:2,3:4) = dt*eye(2);
% B = [zeros(2); dt*eye(2)]; 
B = [dt^2/2 0; 0 dt^2/2; dt 0; 0 dt];
% C = eye(n);
C = [1 0 0 0; 0 1 0 0];
K = dlqr(A,B,eye(n),eye(2)); % feedback law u = -Kx
% Q = 0.001 * eye(n); % process noise covariance
Q = 0.1*[dt^3/3 0  dt^2/2 0; 0 dt^3/3 0 dt^2/2; dt^2/2 0 dt 0; 0 dt^2/2 0 dt];
% R = 0.001 * eye(n); % measurement noise covariance
R = 0.001*eye(2);

% form sys struct
sys.A = A; sys.B = B; sys.C = C;
sys.K = K; sys.Q = Q; sys.R = R;

x0 = zeros(n,1); % inital state
% P0 = 0.003*eye(n); % initial covariance
P0 = diag([0.01 0.01 0.0001 0.0001]);

% define trajectory parameter space
K1 = 1;
K2 = 1; 
U_nom = zonotope([0 K1 0; 0 0 K2]);
% U_nom = zonotope([1 0 0; 1 0 0]);

%% full trajectory reachability calculation

% initial reachable set X0
% xZ0 = [0;0;0;0]; % (uncertain) mean
%                  % x0, y0, k1_c, k2_c
% xZ0 = [xZ0, [0;0;K1;0], [0;0;0;K2]]; % K1, K2 generators
X0 = probZonotope(x0,cov2probGen(P0),3);
%X0 = probZonotope(xZ0,cov2probGen(blkdiag(P0,zeros(2))),3);

tic
pXrs = LQG_matrix_FRS(U_nom,x0,sys,X0,N);
toc

%% plot 3-sigma confidence zonotopes
m = 3;
Xrs = cell(1,N);

figure(1)
hold on; grid on;

% plot full FRS
for i = 1:N
    umeanZ = mean(pXrs{i});
    covZ = cov2zonotope(sigma(pXrs{i}),m,n);
%     covZ = cov2zonotope(sigma(pXrs{i}),m);
    Xrs{i} = umeanZ + covZ;
    plot(Xrs{i},[1,2],'Color','black');
end

axis equal
xlabel('x-coordinate (m)');
ylabel('y-coordinate (m)');

%% create obstacle and intersect with FRS
obs_center = [5; 0];
obs_gen = [[1; 0], [0.5; 1.5]];
obs_zono = zonotope([obs_center, obs_gen]);
obstacle = obs_zono.Z;

% plot obstacle;
figure(1);
p_obs = plot(obs_zono,[1, 2],'FaceColor','r','Filled',true);
p_obs.FaceAlpha = 0.5;
p_obs.EdgeAlpha = 0.5;

% loop through zonotopes of FRS, find "k-sliceable" generators
% and generate obstacle constraints.
obs_dim = [1; 2]; % note that the obstacle exists in the x-y space (not theta or v)
k_dim = [5; 6]; % note that the parameters k are in the 3rd and 4th rows of the zonotopes
buffer_dist = 0; % assume no buffer.

% polytope (halfspace) representation of obstacle/FRS intersection
% this allows us to formulate constraints for the optimization 
A_con = {};
b_con = {};
for i = 1:length(Xrs)
    Z = Xrs{i}.Z;
    c = Z(obs_dim, 1); 
    G = Z(:, 2:end);
    
    [~, k_col] = find(G(k_dim, :) ~= 0); % find "k-sliceable" generators
    k_slc_G = G(obs_dim, k_col);
    k_no_slc_G = G(obs_dim, :); 
    k_no_slc_G(:, k_col) = [];
    
    buff_obstacle_c = [obstacle(:, 1) - c];
    buff_obstacle_G = [obstacle(:, 2:end), k_no_slc_G, buffer_dist*eye(2)]; % obstacle is "buffered" by non-k-sliceable part of FRS
    buff_obstacle_G(:, ~any(buff_obstacle_G)) = []; % delete zero columns of G
    buff_obstacle = [buff_obstacle_c, buff_obstacle_G];
    [A_obs, b_obs] = polytope_PH(buff_obstacle); % turn zonotope into polytope
    
    A_con{i} = A_obs*k_slc_G; % now polytope is over coefficients of k_slc_G generators
    b_con{i} = b_obs;
end

% plot parameters and unsafe set
figure(2); clf; hold on;
title('Set of parameters (red is unsafe)', 'FontSize', 24);
xlabel('$k_1$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$k_2$', 'Interpreter', 'latex', 'FontSize', 24);

% get vertices of k limits for plotting:
lims = [-K1, -K1, K1, K1, -K1; -K2, K2, K2, -K2, -K2];
plot(lims(1, :)', lims(2, :)', 'k', 'LineWidth', 4);

% grid over parameter space
k1_sample = linspace(-K1, K1, 50);
k2_sample = linspace(-K2, K2, 50);
[Xk, Yk] = meshgrid(k1_sample, k2_sample);
c_k = [0;0]; g_k = [K1;K2];

Zk = inf*ones(size(Xk));
for i = 1:length(k1_sample)
    for j = 1:length(k2_sample)
        K = [Xk(i, j); Yk(i, j)];
        lambdas = (K - c_k)./g_k; % given a parameter, get coefficients on k_slc_G generators
        for k = 2:length(A_con) % ============ since initial set is not sliceable w.r.t. inputs
            Zk_tmp = A_con{k}*lambdas - b_con{k}; % A*lambda - b <= 0 means inside unsafe set
            Zk_tmp = max(Zk_tmp); % max of this <= 0 means inside unsafe set
            Zk(i, j) = min(Zk(i, j), Zk_tmp); % take smallest max. if it's <=0, then unsafe
        end
    end
end
p_unsafe = contourf(Xk, Yk, -Zk, [0, 0], 'FaceColor', 'r'); % show zero level set contours

%% have user select a point, see the corresponding slice of FRS:
figure(2);
disp('Click a point!');
[k1_user, k2_user] = ginput(1);
plot(k1_user, k2_user, 'b.', 'MarkerSize', 30, 'LineWidth', 6);

% plot reachable set corresponding to particular parameter "slice"
figure(1);
% ============ since initial set is not sliceable w.r.t. inputs
% p_0 = plot(X0,[1, 2],'FaceColor','b','Filled',true);
% p_0.FaceAlpha = 0.25;
for i = 2:length(Xrs)
    p_slice = plot(zonotope_slice(Xrs{i},k_dim,[k1_user; k2_user]),[1, 2],'FaceColor','b','Filled',true);
    p_slice.FaceAlpha = 0.25;
end

%% sample trajectories 

N_traj = 1000;

u = repmat([k1_user; k2_user],1,N);
X = simulate_LQG_trajectory(sys,N,N_traj,u,x0,P0);

figure(1);
for i_traj = 1:N_traj
    plot(X(1,:,i_traj),X(2,:,i_traj),'r');
end