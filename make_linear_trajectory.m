function [T,U,Z] = make_linear_trajectory(t_f,vx_des,vy_des)
% [T,U,Z] = make_turtlebot_desired_trajectory(t_f,w_des,v_des)
%
% Create a Dubins path as a full-state trajectory for the TurtleBot.
%
% The inputs are:
%   t_f       planning time horizon
%   vx_des    desired x speed
%   vy_des    desired y speed
%
% The outputs are:
%   T        timing for desired trajectory as a 1-by-N array
%   U        desired input (x and y speed) as 2-by-N array
%   Z        desired trajectory (x,y,vx,vy) as a 4-by-N array
%
% Author: Adam Dai
% Created: 14 July 2020

    % set up timing
    t_sample = 0.01 ;
    T = unique([0:t_sample:t_f,t_f],'stable');
    N_t = length(T) ;
    
    % get inputs for desired trajectories
    vx_traj = vx_des*ones(1,N_t) ;
    vy_traj = vy_des*ones(1,N_t) ;
    U_in = [vx_traj ; vy_traj] ;

    % compute desired trajectory
    z0 = zeros(2,1) ;
    [~,Z] = ode45(@(t,z) linear_trajectory_producing_model(t,z,T,U_in),T,z0) ;
    Z = Z';
    
    % compute inputs for robot
    U = U_in ;
end

function zd = linear_trajectory_producing_model(t,z,T_in,U_in)
% zd = linear_trajectory_producing_model(t,z,T_in,U_in)
%
% Output the dynamics of linear 2D point motion
%
% Author: Adam Dai
% Created: 14 July 2020

    % get inputs
    u = match_trajectories(t,T_in,U_in) ;
    vx_des = u(1) ;
    vy_des = u(2) ;

    % compute dynamics
    zd = [vx_des; vy_des] ;
end