function zd = linear_symbolic_traj_prod_model(z,vx,vy)
% zd = turtlebot_symbolic_traj_prod_model(z,v_des,w_des)
%
% Output the dynamics of the linear agent's trajectory-producing model
% that can be used with symbolic variables.
%

%     % extract states
%     x = z(1) ;
%     y = z(2) ;
    
    % compute dynamics
    zd = [vx ; vy ] ;
end