% example for computing FRS for linear agent's (nominal) trajectory
% and intersecting it with a zonotopic obstacle

%% compute reachable set
clear; clc;

t_total = 1;
dt = 0.01;

% model initial condition set:
% X0 = [x_0; % initial x position
%       y_0; % initial y position
%       K1;  % desired x speed (K1 parameter set, centered at k1c with generator k1g)
%       K2]; % desired y speed (K2 parameter set, centered at k2c with generator k2g)

% non-sliceable initial conditions
x_0 = 0;
y_0 = 0;

% sliceable parameters... desired x speed and y speed
k1c = 0; % center of x speed interval
k1g = 0.2; % k1c +- k1g
k2c = 0; % center of y speed interval
k2g = 0.2; % k2c +- k2g
c_k = [k1c; k2c];
g_k = [k1g; k2g];

% set FRS options
options.tStart = 0;
options.tFinal = t_total;
options.x0 = [x_0; y_0; k1c; k2c]; % center of initial set
options.R0 = zonotope([options.x0, [0;0;k1g;0], [0;0;0;k2g]]); % generators for parameter dimensions
options.timeStep=dt; %time step size for reachable set computation
options.taylorTerms=5; %number of taylor terms for reachable sets
options.zonotopeOrder= 20; %zonotope order... increase this for more complicated systems.
options.maxError = 1*ones(4, 1); % this controls splitting, set it high to avoid splitting
options.verbose = 1;

options.uTrans = 0; % center of input set
options.U = zonotope([0, 0]);

options.originContained = 1;
options.advancedLinErrorComp = 0;
options.tensorOrder = 1;
options.reductionInterval = inf;
options.reductionTechnique = 'girard';

% specify continuous dynamics----------------------------------------------
% transition matrix for state augmented with constant input (k1,k2)
A = [0 0 1 0; 
     0 0 0 1; 
     0 0 0 0; 
     0 0 0 0]; 
B = zeros(4,1);
sys = linearSys('sys',A,B);

% compute reachable set----------------------------------------------------
tic
Rcont = reach(sys, options);
tComp = toc;
disp(['computation time of reachable set: ', num2str(tComp)]);

% plot reachable set
figure(1); clf; hold on; axis equal;
title('Workspace', 'FontSize', 24);
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 24);
N = length(Rcont);
for i = 1:N
   Rcont{i} = deleteZeros(Rcont{i});
   p_FRS = plotFilled(Rcont{i}, [1, 2], 'g'); 
   p_FRS.FaceAlpha = 0.02;
   p_FRS.EdgeAlpha = 0.4;
end

%% create obstacle and intersect with FRS
obs_center = [0.1; 0];
obs_gen = [[0.03; 0], [0.05; 0.05]];
obs_zono = zonotope([obs_center, obs_gen]);
obstacle = obs_zono.Z;

% plot obstacle;
figure(1);
p_obs = plotFilled(obs_zono, [1, 2], 'r');
p_obs.FaceAlpha = 0.5;
p_obs.EdgeAlpha = 0.5;

% loop through zonotopes of FRS, find "k-sliceable" generators
% and generate obstacle constraints.
obs_dim = [1; 2]; % note that the obstacle exists in the x-y space (not theta or v)
k_dim = [3; 4]; % note that the parameters k are in the 3rd and 4th rows of the zonotopes
buffer_dist = 0; % assume no buffer.

% polytope representation of obstacle
A_con = cell(1,N);
b_con = cell(1,N);
for i = 1:N
    Z = Rcont{i}.Z;
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
lims = [k1c - k1g, k1c - k1g, k1c + k1g, k1c + k1g, k1c - k1g; k2c - k2g, k2c + k2g, k2c + k2g, k2c - k2g, k2c - k2g];
plot(lims(1, :)', lims(2, :)', 'k', 'LineWidth', 4);

% grid over parameter space
k1_sample = linspace(k1c - k1g, k1c + k1g, 50);
k2_sample = linspace(k2c - k2g, k2c + k2g, 50);
[Xk, Yk] = meshgrid(k1_sample, k2_sample);

Zk = inf*ones(size(Xk));
for i = 1:length(k1_sample)
    for j = 1:length(k2_sample)
        K = [Xk(i, j); Yk(i, j)];
        lambdas = (K - c_k)./g_k; % given a parameter, get coefficients on k_slc_G generators
        for k = 1:length(A_con)
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
for i = 1:N
    p_slice = plotFilled(zonotope_slice(Rcont{i}, k_dim, [k1_user; k2_user]), [1, 2], 'b');
    p_slice.FaceAlpha = 0.25;
end
