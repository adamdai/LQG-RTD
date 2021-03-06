% Compute forward stochastic reachable set of an n-dimensional LQG system 
% for a space of trajectory paramters
% 
% Inputs:
%   S     - zonotope representing trajectory parameter space
%           contains v_0_max, a_0_max, and v_peak_max for each x and y
%           dimension
%   sys   - struct containing system matrices A,B,and C, LQR gain K, and  
%           noise covariance matrices Q and R  
%   X0    - initial stochastic set stored as a probZonotope
%   N     - trajectory length
%   LPM   - linear planning model
%
% Outputs:
%   pXrs - cell array of probZonotopes representing forward reachable sets
%          in time
%   

function pXrs = LQG_param_FRS(S, x_nom_0, sys, X0, N, M)
    
    % retrieve system matrices from sys struct
    A = sys.A; B = sys.B; C = sys.C;
    K = sys.K; Q = sys.Q; R = sys.R;
    
    % parameters
    n = size(A,1); % system dimension
    m = dim(S); % trajectory parameter dimension
    
    % initial state covariance
    P = sigma(X0); 
    P = P(1:n,1:n);
    
    % initial state estimation error
    WpZ = probZonotope(zeros(n,1),cov2probGen(Q),3);
    VpZ = probZonotope(zeros(2,1),cov2probGen(R),3);

    % initialize FRS
    pXrs = cell(1,N);
    pXrs{1} = [eye(n); zeros(m,n)]*x_nom_0 + [zeros(n,m); zeros(m,m)] * S + [eye(n); zeros(m,n)]*(X0 + -x_nom_0);

    % recursive coefficients
    coeff_a = eye(n); coeff_b = 0; % initialize 
    coeff_c = cell(1,N); coeff_d = cell(1,N);
    coeff_c{1} = nan; coeff_d{1} = nan;

    coeff_e = eye(n);
    coeff_p = cell(1,N); coeff_q = cell(1,N);
    coeff_p{1} = nan; coeff_q{1} = nan;
    
    
    for t = 2:N
        
        % update coeffs a and b
        coeff_a = (A-B*K)*coeff_a;
        coeff_b = (A-B*K)*coeff_b - B*K*coeff_e;

        % update coeffs c and d
        for i = 2:t-1
            coeff_c{i} = (A-B*K)*coeff_c{i} - B*K*coeff_p{i};
            coeff_d{i} = (A-B*K)*coeff_d{i} - B*K*coeff_q{i};
        end

        % add new coeffs c and d
        coeff_c{t} = eye(n);  coeff_d{t} = zeros(n,2);

        % calculate all CpZ and DpZ terms
        all_CpZ = [coeff_c{t}; zeros(m,n)] * WpZ;
        all_DpZ = [coeff_d{t}; zeros(m,2)] * VpZ;
        for i = 2:t-1
            all_CpZ = all_CpZ + [coeff_c{i}; zeros(m,n)] * WpZ;
            all_DpZ = all_DpZ + [coeff_d{i}; zeros(m,2)] * VpZ;
        end

        % compute reachable set
        AB_coeff = 0;
        for i = 1:t-1
            AB_coeff = AB_coeff + A^(t-i) * B * M{i};
        end

        pXrs{t} = [A^t; zeros(m,n)]*x_nom_0 + [AB_coeff; eye(m,m)] * S + [coeff_a-coeff_b; zeros(m,n)]*(X0 + -x_nom_0) + all_CpZ + all_DpZ;
        
        % online filter steps
        P_pred = A*P*A' + Q;
        L = P_pred*C'/(C*P_pred*C' + R);
        P = P_pred - L*C*P_pred;

        % update coeff e
        coeff_e = (eye(n) - L*C)*A*coeff_e;

        % update coeffs p and q
        for i = 2:t-1
            coeff_p{i} = (eye(n) - L*C)*A*coeff_p{i};
            coeff_q{i} = (eye(n) - L*C)*A*coeff_q{i};
        end

        % add coeffs p and q for new w and v
        coeff_p{t} = -(eye(n) - L*C);  coeff_q{t} = L;

    end
end