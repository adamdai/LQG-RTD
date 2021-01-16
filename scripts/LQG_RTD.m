% Perform RTD over a single planning horizon for an LQG system initialized
% near a single obstacle

clc
clear
close all

%% double integrator system
n_s = 4; % state dimension
n_i = 2; % input dimension
n_m = 2; % measurement dimension
pos_i = 1:2; % position indices

% trajectory discretization and length (planning horizon)
t_f = 3; dt = 0.1; N = t_f/dt + 1; t_peak = 1.5;

% system
A = eye(n_s);
A(1:2,3:4) = dt*eye(2);
B = [dt^2/2 0; 
     0      dt^2/2; 
     dt     0; 
     0      dt]; 
C = [1 0 0 0;  
     0 1 0 0]; 
 
K = dlqr(A,B,eye(n_s),eye(n_i)); % feedback law u = -Kx

Q = 0.01*[dt^3/3 0      dt^2/2 0;  
         0      dt^3/3 0      dt^2/2; 
         dt^2/2 0      dt     0; 
         0      dt^2/2 0      dt]; % process noise covariance
R = 0.0001*eye(n_m); % measurement noise covariance

% form sys struct
sys.A = A; sys.B = B; sys.C = C;
sys.K = K; sys.Q = Q; sys.R = R;

% initial conditions
v0_x = 2;
v0_y = 0;
a0_x = 2;
a0_y = 0;

x0 = [5;5;v0_x;v0_y]; % inital state
P0 = 0.1*diag([0.01 0.01 0.0001 0.0001]); % initial covariance

% define trajectory parameter space
v_peak_max = 3;
S = zonotope([v0_x 0          0; 
              a0_x 0          0; 
              0    v_peak_max 0;
              v0_y 0          0; 
              a0_y 0          0; 
              0    0          v_peak_max]);
          
% construct mapping matrix from LPM
load quadrotor_linear_planning_model.mat
a_mat = LPM.acceleration;
M = {};
for i = 1:N
    M{i} = [a_mat(:,i)' zeros(1,3); 
            zeros(1,3)  a_mat(:,i)'];
end

%% plot agent

figure(1)
hold on; grid on;

scatter(x0(1), x0(2),'filled');

axis equal
xlabel('x-coordinate (m)');
ylabel('y-coordinate (m)');

%% compute FRS

X0 = probZonotope(x0,cov2probGen(P0),3);

tic
pXrs = LQG_param_FRS(S,x0,sys,X0,N,M);
disp('FRS computation:'); toc

%% plot 3-sigma confidence zonotopes
m = 3;
Xrs = cell(1,N);

figure(1)
hold on; grid on;

% plot full FRS
for i = 1:N
    umeanZ = mean(pXrs{i});
    covZ = cov2zonotope(sigma(pXrs{i}),m,n_s);
    Xrs{i} = deleteZeros(umeanZ + covZ);
    plot(Xrs{i},[1,2],'Color','black');
end

axis equal
xlabel('x-coordinate (m)');
ylabel('y-coordinate (m)');


%% create obstacles and generate constraints
obs = {};

obs_center = [3; -3];
obs_gen = [[1; 0], [0; 1]];
obs{1} = zonotope([obs_center, obs_gen]);

obs_center = [2; 3];
obs_gen = [[1; 0], [0; 1]];
obs{2} = zonotope([obs_center, obs_gen]);

% plot obstacles
figure(1);
for i = 1:length(obs)
    p_obs = plot(obs{i},[1, 2],'FaceColor','r','Filled',true);
    p_obs.FaceAlpha = 0.5;
    p_obs.EdgeAlpha = 0.5;
end

% loop through zonotopes of FRS, find "k-sliceable" generators
% and generate obstacle constraints.
obs_dim = [1; 2]; % note that the obstacle exists in the x-y space (not theta or v)
k_dim = [7; 10]; % note that the parameters k are in the 3rd and 4th rows of the zonotopes

% polytope (halfspace) representation of obstacle/FRS intersection
% this allows us to formulate constraints for the optimization 
[A_con, b_con] = generate_obs_constraints(Xrs, obs, k_dim, obs_dim);

%% trajectory optimization
x_g = [6;-2]; % arbitary goal
ic = [v0_x; v0_y; a0_x; a0_y]; % initial conditions
c_k = [0;0]; g_k = [v_peak_max;v_peak_max];
cost = @(k) traj_opt_cost(sys,pos_i,N,k,ic,x_g,M); % cost function
cons = @(k) traj_opt_cons(k,c_k,g_k,A_con,b_con); % nonlinear constraint function 
%cons = [];
initial_guess = zeros(2,1);

lb = [-v_peak_max, -v_peak_max];
ub = [v_peak_max, v_peak_max];

o = optimoptions('fmincon');
o.OptimalityTolerance = 1e-8;
%o.ConstraintTolerance = 1e-10;
o.MaxIterations = 100000;
o.MaxFunctionEvaluations = 10000;
o.SpecifyConstraintGradient = true;
o.SpecifyObjectiveGradient = true;

% fmincon solution
tic
k = fmincon(cost,initial_guess,[],[],[],[],lb,ub,cons,o);
disp('fmincon:'); toc
disp(k)

% unconstrained least-squares solution
% k = ls_k_solver(sys,pos_i,N,ic,x_g,M,lb,ub);
% disp(k)

%% plot workspace (trajectories and reachable sets)

% plot reachable set corresponding to particular parameter "slice"
figure(1);
% ============ since initial set is not sliceable w.r.t. inputs
p_init = plot(Xrs{1},[1, 2],'FaceColor','b','Filled',true);
p_init.FaceAlpha = 0.25;
for i = 2:length(Xrs)
    p_slice{i} = plot(zonotope_slice(Xrs{i},k_dim,k),[1, 2],'FaceColor','b','Filled',true);
    p_slice{i}.FaceAlpha = 0.25;
end



% show goal
scatter(x_g(1), x_g(2),'filled');

% construct control input from traj param and LPM
s = [v0_x; a0_x; k(1); v0_y; a0_y; k(2)];
u = zeros(2,N);
for i = 1:N
    u(:,i) = M{i}*s;
end

% simulate rollouts
N_traj = 100;
[X,x_nom] = simulate_LQG_trajectory(sys,N,N_traj,u,x0,P0);

figure(1);

% plot nominal trajectory
plot(x_nom(1,:),x_nom(2,:),'g');

% plot sampled trajectories
for i_traj = 1:N_traj
    plot(X(1,:,i_traj),X(2,:,i_traj),'r');
end

%% plot parameter space (unsafe set and optimal choice of k)
figure(2); clf; hold on;
title('Set of parameters (red is unsafe)', 'FontSize', 24);
xlabel('$k_1$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$k_2$', 'Interpreter', 'latex', 'FontSize', 24);

% get vertices of k limits for plotting:
lims = [-v_peak_max, -v_peak_max, v_peak_max, v_peak_max, -v_peak_max; 
        -v_peak_max, v_peak_max, v_peak_max, -v_peak_max, -v_peak_max];
plot(lims(1, :)', lims(2, :)', 'k', 'LineWidth', 4);

% grid over parameter space
k1_sample = linspace(-v_peak_max, v_peak_max, 50);
k2_sample = linspace(-v_peak_max, v_peak_max, 50);
[Xk, Yk] = meshgrid(k1_sample, k2_sample);

% determine unsafe set by evaluating sampled k's on constraints
Zk = inf*ones(size(Xk));
Zc = inf*ones(size(Xk));
for i = 1:length(k1_sample)
    for j = 1:length(k2_sample)
        K = [Xk(i, j); Yk(i, j)];
        Zk(i,j) = cons(K);
        Zc(i,j) = cost(K);
    end
end

% plot unsafe set
p_unsafe = contourf(Xk, Yk, Zk, [0, 0], 'FaceColor', 'r'); % show zero level set contours
%cost_contours = contourf(Xk, Yk, Zc, 100); % show cost function contours

% show optimal choice of k
scatter(k(1), k(2),'filled','b');
