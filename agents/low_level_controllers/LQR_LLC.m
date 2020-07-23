classdef LQR_LLC < low_level_controller
    properties

        interp_mode = 'previous' ;
        discrete = 0;  % discrete or continuous time LQR
        estimating = 0;  % whether to apply feedback based on true state
                         % or state estimate
        Q = eye(2);
        R = eye(2);
        
    end
    
    methods
        %% constructor
        function LLC = LQR_LLC(varargin)
            n_agent_states = 2 ;
            n_agent_inputs = 2 ;
            
            LLC = parse_args(LLC,'n_agent_states',n_agent_states,...
                'n_agent_inputs',n_agent_inputs,varargin{:}) ;
        end
        
        %% get control inputs
        function U = get_control_inputs(LLC,A,t_cur,z_cur,T_des,U_des,Z_des)
            % get current state
            if LLC.estimating
                p_cur = z_cur(A.state_est_indices) ;
            else
                p_cur = z_cur(A.position_indices) ;
            end
            
            % get desired state and inputs
            [u_des,z_des] = match_trajectories(t_cur,T_des,U_des,T_des,Z_des,LLC.interp_mode) ;
            p_des = z_des(A.position_indices) ;
            
            % LQR feedback gain K
            if LLC.discrete
                K = dlqr(A.A,A.B,LLC.Q,LLC.R);
            else
                K = lqr(A.A,A.B,LLC.Q,LLC.R);
            end
            U = u_des - K*(p_cur - p_des);
        end
    end
end