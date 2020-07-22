classdef discrete_agent < translating_agent_2D
% Class: linear_agent < translating_agent_2D < agent
%
% Description
    
    properties        
        % system model
        A = eye(2); 
        B = 0.01*eye(2);     
        C = eye(2);
        Q = 0.001*eye(2);
        R = 0.001*eye(2);
        
        % state estimate indices
        state_est_indices = [3 4];
        cov_indices = [5 6 7 8];
        
        % input limits 
        max_xspeed = 2.0 ; % m/s  
        max_yspeed = 2.0 ; % m/s 
        
        % integrator type, to allow for fixed time step integration
        integrator_type = 'discrete' ; % choose 'ode45' or 'ode4' or 'ode113'
        integrator_time_discretization = 0.01 ; % for ode4
    end
    
    methods
        %% constructor
        function A = discrete_agent(varargin)
            % set up default superclass values
            name = 'Discrete Agent' ;
            default_footprint = 0.35/2 ;
            n_states = 8 ;
            n_inputs = 2 ;
            stopping_time = 3 ; % conservative estimate
            sensor_radius = 3 ;
            LLC = LQR_LLC ;
            
            % create agent
            A@translating_agent_2D('name',name,...
                'footprint',default_footprint,...
                'n_states',n_states,'n_inputs',n_inputs,...
                'stopping_time',stopping_time,'sensor_radius',sensor_radius,...
                'LLC',LLC,varargin{:}) ;
        end
        
        %% emergency stop
        % note, this ignores any previous trajectory the agent may have
        % been tracking; we have to define this different from the default
        % for the TurtleBot because it takes in acceleration as a control
        % input, as opposed to doing feedback about desired speed
        function stop(A,t_stop)
            if nargin < 2
                t_stop = A.stopping_time ;
            end
            
            % get the current speed
            v = A.state(A.speed_index,end) ;

            % check how long it will take to come to a stop and make the
            % stopping time vector
            t_req_to_stop = v/A.max_accel ;            
            T_input = [0, max(t_req_to_stop,t_stop)] ;
            
            % generate the input and desired trajectory
            U_input = zeros(2,2) ;
            Z_input = [repmat(A.state(1:3,end),1,2); zeros(1,2)] ;
            
            % call move method to perform stop
            % A.LLC.gains = A.LLC.stop_gains ;
            A.move(t_stop,T_input,U_input,Z_input) ;
            
            % reset the default gains after stopping
            % A.LLC.gains = A.LLC.default_gains ;
        end
        
        %% dynamics
        function sn = dynamics(A,t,z,T,U,Z)
            % get current state
            x = z(A.position_indices) ;
            x_est = z(A.state_est_indices) ;
            P = reshape(z(A.cov_indices),2,2) ;
            
            % get feedback control inputs
            u = A.LLC.get_control_inputs(A,t,z,T,U,Z) ;
            
            % saturate the inputs
            u(1) = bound_values(u(1),A.max_xspeed) ;
            u(2) = bound_values(u(2),A.max_yspeed) ;
            
            % sample noise
            w = mvnrnd([0;0],A.Q)';
            v = mvnrnd([0;0],A.R)';
            
            % update state
            xn = A.A*x + A.B*u + w;
            y = A.C*xn + v;
            
            x_estn = A.A*x_est + A.B*u;
            Pn = A.A*P*A.A' + A.Q;
            
            L = Pn*A.C'/(A.C*Pn*A.C' + A.R);
            x_estn = x_estn + L*(y - A.C*x_estn);
            Pn = Pn - L*A.C*Pn;
            sn = [xn; x_estn; vec(Pn)];
        end
        
        %% integrator options
        function [tout,zout] = integrator(A,fun,tspan,z0)
            switch A.integrator_type
                case 'ode45'
                    [tout,zout] = ode45(@(t,z) fun(t,z),tspan,z0(:)) ;
                case 'ode113'
                    [tout,zout] = ode113(@(t,z) fun(t,z),tspan,z0(:)) ;
                case {'ode4','RK4'}
                    dt = A.integrator_time_discretization ;
                    tout = tspan(1):dt:tspan(end) ;
                    if tout(end) ~= tspan(end)
                        tout = [tout, tspan(end)] ;
                    end
                    zout = ode4(@(t,z) fun(t,z),tout,z0(:)) ;
                case 'discrete'
                    dt = A.integrator_time_discretization ;
                    tout = tspan(1):dt:tspan(end) ;
                    if tout(end) ~= tspan(end)
                        tout = [tout, tspan(end)] ;
                    end
                    N = length(tout);
                    zout = zeros(A.n_states,N);
                    zout(:,1) = z0;
                    for i = 2:N
                        zout(:,i) = fun(tout(i-1),zout(:,i-1));
                    end
                    zout = zout';
                otherwise
                    error('Please set A.integrator_type to either ode45 or ode4')
            end
            tout = tout(:)' ;
            zout = zout' ;
        end
    end
end