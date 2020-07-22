classdef LQG_LLC < low_level_controller
    properties

        interp_mode = 'previous' ;
        Q = eye(2);
        R = eye(2);
        
    end
    
    methods
        %% constructor
        function LLC = linear_LLC(varargin)
            n_agent_states = 2 ;
            n_agent_inputs = 2 ;
            
            LLC = parse_args(LLC,'n_agent_states',n_agent_states,...
                'n_agent_inputs',n_agent_inputs,varargin{:}) ;
        end
        
        %% get control inputs
        function U = get_control_inputs(LLC,A,t_cur,z_cur,T_des,U_des,Z_des)
            % get current state estimate
            p_est = z_cur(A.state_est_indices) ;
            
            % get desired state and inputs
            [u_des,z_des] = match_trajectories(t_cur,T_des,U_des,T_des,Z_des,LLC.interp_mode) ;
            p_des = z_des(A.position_indices) ;
            
            % LQR feedback gain K
            K = lqr(A.A,A.B,LLC.Q,LLC.R);
            U = u_des - K*(p_est - p_des);
        end
    end
end