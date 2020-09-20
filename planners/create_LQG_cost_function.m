%% description
% In this script, we create a cost function and its analytic derivative for
% the LQG agent. The cost function is simply to minimize distance to a
% waypoint while tracking a trajectory. We'll save the cost and gradient as
% MATLAB functions.
%
%% user parameters
% numerical integration parameters
dt_int = 0.01 ; % s
discount_factor = [1.0 ; 1.0] ; % compensates for integration error in x/y

% % desired params for validation
k_1_val = 0.5 ;
k_2_val = 0.5 ;

% save cost and gradient flag
save_flag = true ;

%% automated from here
% load timing
load('turtlebot_timing.mat')

% create symbolic variables to use in the cost
syms k_1 k_2 real

%% numerically compute endpoint of desired traj
disp('Numerically estimating desired traj. endpoint')

% we'll express the endpoint of the desired trajectory at t_plan using
% numerical integration - first, set up the initial condition
z = [0;0] ; % (x,y)

tmax = 300 ; % timer just in case
tcur = tic ;
for tidx = 0:dt_int:t_plan
    dzdt = linear_symbolic_traj_prod_model(z,k_1,k_2) ;
    z = z + discount_factor.*dt_int.*dzdt ;

    if toc(tcur) > tmax
        break
    end
end
toc(tcur)

% make a function out of the endpoint to use for validation
z_fn = matlabFunction(z,'Vars',[k_1,k_2]) ;

%% make the cost and analytic gradient
% create waypoints
syms x_des y_des real

% make the cost
cost = sum(([x_des;y_des] - z).^2) ;

% get the analytic gradient
cost_grad = [diff(cost,k_1), diff(cost,k_2)] ;

%% make a function out of the cost and gradient
if save_flag
    disp('Saving cost and cost gradient')
    matlabFunction(cost,'Vars',[k_1,k_2,x_des,y_des],'File','LQG_cost')
    matlabFunction(cost_grad,'Vars',[k_1,k_2,x_des,y_des],'File','LQG_cost_grad')
else
    disp('Not saving cost and cost gradient')
end

%% time the cost and gradient
try
    test_input = num2cell(2*rand(1,6)) ;
    disp('Timing cost function')
    cost_time = timeit(@() turtlebot_cost(test_input{:})) ;
    disp(['    Average exec time: ', num2str(cost_time),' s'])
    
    disp('Timing gradient function')
    cost_grad_time = timeit(@() turtlebot_cost_grad(test_input{:})) ;
    disp(['    Average exec time: ', num2str(cost_grad_time),' s'])
catch
    disp('Cost and gradient not found on path!')
end

%% validation
% get the actual endpoint at t_plan
[~,~,Z] = make_linear_trajectory(t_plan,k_1_val,k_2_val) ;

z_true = Z(1:2,end) ;

% sub in k_1, and k_2 to z
z_val = double(subs(z,[k_1,k_2], [k_1_val,k_2_val])) ;

%% plotting
figure(1) ; clf ; hold on ; axis equal ; grid on

% plot desired trajectory
plot(Z(1,:),Z(2,:),'b')

% plot approximate endpoint
plot(z_val(1),z_val(2),'bo') ;

% cleanup
legend('desired traj.','numerical endpoint')

set(gca,'FontSize',15)
xlabel('x')
ylabel('y')