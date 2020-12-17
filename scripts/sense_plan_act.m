%% SPA simulation (Sense-Plan-Act)

% time intervals
t_integration = 0.01 ;
t_IMU = 0.01 ;
t_LIDAR = 0.1 ;
t_plan = 1.0 ;
t_controller = 0.1 ; 


tspan = [0 4] ;
x0 = [0; 0; 0] ;
[t,x] = ode45(@(t,x) f, tspan, x0) ; 

% dynamics: Dubin's
% x(0) - x      u(0) - V
% x(1) - y      u(1) - omega
% x(2) - theta
function x_dot = dynamics(x,u)
    x_dot(0) = u(0) * cos(x(2)) ; 
    x_dot(1) = u(0) * sin(x(2)) ; 
    x_dot(2) = u(1) ; 
end

% controller
function u = controller(t,x)

end

