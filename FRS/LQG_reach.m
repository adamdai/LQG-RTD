% Compute forward stochastic reachable set of an n-dimensional LQG system 
% along a nominal trajectory of length N given an initial stochastic set
% 
% Inputs:
%   x_nom - n x N array containing nominal trajectory, or cell array of 
%           zonotopes containing nominal trajectory reachable space
%   sys   - struct continaing system matrices A,B,and C, LQR gain K, and  
%           noise covariance matrices Q and R  
%   X0    - initial stochastic set stored as a probZonotope
%
% Outputs:
%   pXrs - cell array of probZonotopes representing forward reachable sets
%          in time
%   

function pXrs = LQG_reach(x_nom, sys, X0)
    
    % retrieve system matrices from sys struct
    A = sys.A; B = sys.B; C = sys.C;
    K = sys.K; Q = sys.Q; R = sys.R;
    
    N = size(x_nom,2);
    P = sigma(X0);
    
    % process noise and measurement noise
    WpZ = probZonotope([0;0],cov2probGen(Q),3);
    VpZ = probZonotope([0;0],cov2probGen(R),3);

    % initialize FRS
    pXrs = cell(1,N);
    pXrs{1} = X0;
    
    % initial state estimation error set
    X_err = X0 + -center(X0);

    % coeff_a = (A-B*K); coeff_b = -B*K; 
    coeff_a = eye(2); coeff_b = 0; % initialize 
    coeff_c = cell(1,N); coeff_d = cell(1,N);
    coeff_c{1} = nan; coeff_d{1} = nan;

    coeff_e = eye(2);
    coeff_p = cell(1,N); coeff_q = cell(1,N);
    coeff_p{1} = nan; coeff_q{1} = nan;

    for k = 2:N

        % update coeffs a and b
        coeff_a = (A-B*K)*coeff_a;
        coeff_b = (A-B*K)*coeff_b - B*K*coeff_e;

        % update coeffs c and d
        for n = 2:k-1
            coeff_c{n} = (A-B*K)*coeff_c{n} + -B*K*coeff_p{n};
            coeff_d{n} = (A-B*K)*coeff_d{n} + -B*K*coeff_q{n};
        end

        % add new coeffs c and d
        coeff_c{k} = eye(2);  coeff_d{k} = zeros(2);

        % calculate all CpZ and DpZ terms
        all_CpZ = coeff_c{k}*WpZ;
        all_DpZ = coeff_d{k}*VpZ;
        for n = 2:k-1
            all_CpZ = all_CpZ + coeff_c{n}*WpZ;
            all_DpZ = all_DpZ + coeff_d{n}*VpZ;
        end
        
        % online filter steps
        P_pred = A*P*A' + Q;
        L = P_pred*C'/(C*P_pred*C' + R);
        P = P_pred - L*C*P_pred;

        % update coeff e
        coeff_e = (eye(2) - L*C)*A*coeff_e;

        % update coeffs p and q
        for n = 2:k-1
            coeff_p{n} = (eye(2) - L*C)*A*coeff_p{n};
            coeff_q{n} = (eye(2) - L*C)*A*coeff_q{n};
        end

        % add coeffs p and q for new w and v
        coeff_p{k} = -(eye(2) - L*C);  coeff_q{k} = L;
        
        % compute state estimator error
        ILC = (eye(2) - L*C);
        X_err = ILC*(A*X_err + WpZ) + L*VpZ;
        
        % compute reachable set
        %pXr = (coeff_a - coeff_b)*(X0 + -x_nom(:,1)) + all_CpZ + all_DpZ; % assumes initial state estimate equals initial nominal point
        pXr = coeff_a*(X0 + -x_nom(:,1)) + coeff_b*X_err + all_CpZ + all_DpZ; % assumes initial state and initial state estimation error are independent
        pXrs{k} = pXr + x_nom(:,k);

    end
end