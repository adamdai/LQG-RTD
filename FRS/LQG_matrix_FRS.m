% Compute forward stochastic reachable set of an n-dimensional LQG system 
% for a space of trajectory paramters
% 
% Inputs:
%   U_nom - trajectory parameter space represented by a cell array of 
%           zonotopes
%   sys   - struct containing system matrices A,B,and C, LQR gain K, and  
%           noise covariance matrices Q and R  
%   X0    - initial stochastic set stored as a probZonotope
%   N     - trajectory length
%
% Outputs:
%   pXrs - cell array of probZonotopes representing forward reachable sets
%          in time
%   

function pXrs = LQG_matrix_FRS(U_nom, x_nom_0, sys, X0, N)
    
    % retrieve system matrices from sys struct
    A = sys.A; B = sys.B; C = sys.C;
    K = sys.K; Q = sys.Q; R = sys.R;
    
    % parameters
    n = size(A,1); % system dimension
    
    % initial state covariance
    P = sigma(X0); 
    P = P(1:n,1:n);
    
    % initial state estimation error
    X_0_tilde = probZonotope([0;0], cov2probGen(P),3);
    
    % process noise and measurement noise
%     WpZ = probZonotope([0;0;0;0],blkdiag(cov2probGen(Q),zeros(2)),3);
%     VpZ = probZonotope([0;0;0;0],blkdiag(cov2probGen(R),zeros(2)),3);
    WpZ = probZonotope([0;0],cov2probGen(Q),3);
    VpZ = probZonotope([0;0],cov2probGen(R),3);

    % initialize FRS
    pXrs = cell(1,N);
    %pXrs{1} = X0;

    % recursive coefficients
    % coeff_a = (A-B*K); coeff_b = -B*K; 
    coeff_a = eye(n); coeff_b = 0; % initialize 
    coeff_c = cell(1,N); coeff_d = cell(1,N);
    coeff_c{1} = nan; coeff_d{1} = nan;

    coeff_e = eye(n);
    coeff_p = cell(1,N); coeff_q = cell(1,N);
    coeff_p{1} = nan; coeff_q{1} = nan;

    for k = 1:N
        
        % update coeffs a and b
        coeff_a = (A-B*K)*coeff_a;
        coeff_b = (A-B*K)*coeff_b - B*K*coeff_e;

        % update coeffs c and d
        for i = 2:k-1
            coeff_c{i} = (A-B*K)*coeff_c{i} + -B*K*coeff_p{i};
            coeff_d{i} = (A-B*K)*coeff_d{i} + -B*K*coeff_q{i};
        end

        % add new coeffs c and d
        coeff_c{k} = eye(n);  coeff_d{k} = zeros(n);

        % calculate all CpZ and DpZ terms
        all_CpZ = [coeff_c{k}; zeros(2)] * WpZ;
        all_DpZ = [coeff_d{k}; zeros(2)] * VpZ;
        for i = 2:k-1
            all_CpZ = all_CpZ + [coeff_c{i}; zeros(2)] * WpZ;
            all_DpZ = all_DpZ + [coeff_d{i}; zeros(2)] * VpZ;
        end

        % compute reachable set
        AB_coeff = 0;
        for i = 0:k-1
            AB_coeff = AB_coeff + A^(k-1-i);
        end
        AB_coeff = AB_coeff * B;
        
        pXr = [A^k; zeros(2)] * X0 + [AB_coeff; eye(2)] * U_nom + [coeff_a; zeros(2)] * (X0) + [coeff_b; zeros(2)] * X_0_tilde + all_CpZ + all_DpZ;
        %pXr = blkdiag(coeff_a - coeff_b,eye(n))*(X0 - X_nom{1}) + aug_CpZ + aug_DpZ;
        %pXr = blkdiag(coeff_a - coeff_b,eye(n))*X_diff + aug_CpZ + aug_DpZ;
        %pXr = aug_CpZ + aug_DpZ;
        pXrs{k} = pXr;
        
        % online filter steps
        P_pred = A*P*A' + Q;
        L = P_pred*C'/(C*P_pred*C' + R);
        P = P_pred - L*C*P_pred;

        % update coeff e
        coeff_e = (eye(n) - L*C)*A*coeff_e;

        % update coeffs p and q
        for i = 2:k-1
            coeff_p{i} = (eye(n) - L*C)*A*coeff_p{i};
            coeff_q{i} = (eye(n) - L*C)*A*coeff_q{i};
        end

        % add coeffs p and q for new w and v
        coeff_p{k} = -(eye(n) - L*C);  coeff_q{k} = L;

    end
end