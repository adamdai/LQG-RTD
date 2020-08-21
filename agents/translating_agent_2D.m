classdef translating_agent_2D < agent
    %% properties
    properties
        
        % trajectory info
        desired_time
        desired_input
        desired_trajectory
        stopping_time = 1 ;
        
        % footprint info
        footprint = 1 ;
        
        % plotting info
        footprint_vertices
        plot_footprint_color = [0 0 1] ;
        plot_footprint_edge_color = [0 0 1] ;
        plot_footprint_opacity = 0.2 ;
        plot_footprint_edge_opacity = 1 ;
        plot_trajectory_at_time_flag = true ;
    end
    
    %% methods
    methods
        %% constructor
        function A = translating_agent_2D(varargin)
            A@agent('name','Translating Agent 2D',varargin{:}) ;
            
            % set up plot data
            A.plot_data.trajectory = [] ;
            A.plot_data.footprint = [] ;
            
            % set up footprint and arrow for plot
            A.make_footprint_plot_data() ;
            
            % reset time, states, and inputs, just in case
            A.reset() ;
        end
        
        %% reset
        function reset(A,state)
            if nargin < 2
                state = zeros(A.n_states,1) ;
            end
            
            % do the reset
            A.state = zeros(A.n_states,1) ;
            A.time = 0 ;
            A.input = zeros(A.n_inputs,1) ;
            A.input_time = 0 ;
            
            % reset the state
            switch length(state)
                case A.n_states
                    A.state = state ;
                case 2
                    A.state(A.position_indices) = state ;
                otherwise
                    error(['The provided state has an incorrect number of elements!',...
                        ' Please provide either a 2-by-1 position or an n_states-by-1 full state vector.'])
            end
        end
        
        %% get agent into
        function agent_info = get_agent_info(A)
            % call superclass method
            agent_info = get_agent_info@agent(A) ;
            
            % additional fields
            agent_info.desired_time = A.desired_time ;
            agent_info.desired_input = A.desired_input ;
            agent_info.desired_trajectory = A.desired_trajectory ;
            agent_info.footprint = A.footprint ;
            agent_info.footprint_vertices = A.footprint_vertices ;
        end
        
        %% commit move data
        function commit_move_data(A,T_state,Z_state,T_used,U_used,Z_used)
            % commit the usual move data
            commit_move_data@agent(A,T_state,Z_state,T_used,U_used,Z_used) ;
            
            % also update the desired time/input/trajectory
            if ~isempty(A.desired_time)
                t_cur = A.desired_time(end-1) ;
            else
                t_cur = 0 ;
            end
            A.desired_time = [A.desired_time, T_used(2:end) + t_cur, nan(1,1)] ;
            
            % save the desired input
            n_U = size(U_used,1) ;
            A.desired_input = [A.desired_input, U_used(:,2:end), nan(n_U,1)] ;
            
            % save the desired trajectory
            if ~isempty(Z_used)
                n_Z = size(Z_used,1) ;
                A.desired_trajectory = [A.desired_trajectory, Z_used(:,2:end), nan(n_Z,1)] ;
            end
        end
        
        %% emergency stop
        % note, this ignores any previous trajectory the agent may have
        % been tracking
        function stop(A,t_stop)
            if nargin < 2
                t_stop = A.stopping_time ;
            end
            
            T_input = 0:0.1:t_stop ;
            N_t = length(T_input) ;
            
            pose =  A.state([A.position_indices,A.heading_index],end) ;
            stopped_state = [pose ; zeros(A.n_states-3,1)] ;
            
            Z_desired = repmat(stopped_state,1,N_t) ;
            
            U_input = zeros(A.n_inputs,N_t) ;
            
            A.move(t_stop,T_input,U_input,Z_desired) ;
        end
        
        %% utility
        function make_footprint_plot_data(A)
            switch length(A.footprint)
                case 1
                    A.footprint_vertices = make_circle(A.footprint) ;
                case 2
                    A.footprint_vertices = make_box(A.footprint) ;
            end
        end
        
        %% plotting
        function plot(A,~)
            A.plot_at_time(A.time(end)) ;
            plot@agent(A,A.plot_footprint_color) ;
        end
        
        function plot_at_time(A,t)
            % compute footprint for plot
            z_t = match_trajectories(t,A.time,A.state) ;
            p_t = z_t(A.position_indices) ;
            fp_t = A.footprint_vertices(:,1:end-1) ;
            N_fp = size(fp_t,2) ;
            V_fp = fp_t + repmat(p_t,1,N_fp) ;
            
            % plot
            if check_if_plot_is_available(A,'footprint')
                A.plot_data.footprint.Vertices = V_fp' ;
            else
                % plot footprint
                fp_data = patch(V_fp(1,:),V_fp(2,:),A.plot_footprint_color,...
                    'EdgeColor',A.plot_footprint_edge_color,...
                    'FaceAlpha',A.plot_footprint_opacity,...
                    'EdgeAlpha',A.plot_footprint_edge_opacity) ;
                
                % save plot data
                A.plot_data.footprint = fp_data ;
            end
            
            if A.plot_trajectory_at_time_flag
                % get the executed path up to the current time
                X = A.state(A.position_indices,:) ;
                T_log = A.time <= t ;
                X = X(:,T_log) ;
                
                % plot it
                if check_if_plot_is_available(A,'trajectory')
                    A.plot_data.trajectory.XData = X(1,:) ;
                    A.plot_data.trajectory.YData = X(2,:) ;
                end
                    traj_data = plot_path(X,'b-') ;
                    A.plot_data.trajectory = traj_data ;
            end
        end
    end
end