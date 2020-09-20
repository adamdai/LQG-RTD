function [c, gc] = LQG_cost_for_fmincon(k,wp_local,start_tic,timeout)
%
% Evaluate the cost and cost gradient for use with fmincon when performing
% trajectory optimization for the LQG agent.
%

    % evaluate cost and gradient
    c = LQG_cost(k(1),k(2),wp_local(1),wp_local(2)) ;
    gc = LQG_cost_grad(k(1),k(2),wp_local(1),wp_local(2)) ;
    
    % perform timeout check
    if nargin > 3 && toc(start_tic) > timeout
        error('Timed out while evaluating cost function!')
    end
end