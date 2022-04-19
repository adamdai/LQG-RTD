clc
clear
% close all

%% init variables
bhw_x = 0; z_cov_x = 1;
bhw_y = 5;  z_cov_y = 1;

dx = 1; dt = 0.2; N = 50;
A = eye(2); B = dt*eye(2); C = eye(2);
K = dlqr(A,B,eye(2),eye(2));

x0 = [0;0]; P0 = eye(2);
Q = diag([0.1 0.1]);
x_nom = zeros(2,N); x_nom(1,:) = 0:dx:(N-1)*dx;
u_nom = [(dx/dt);0];

N_traj = 100;
p_bias_jump = 0;

R_ob_x = ( ( bhw_x + 3*sqrt(z_cov_x) )/3 )^2;
R_ob_y = ( ( bhw_y + 3*sqrt(z_cov_y) )/3 )^2;
R = diag([R_ob_x R_ob_y]);
mm_available = true;
pZ_order = 1;

%% generate sample trajectories
X = nan(2,N,N_traj);
Z = nan(2,N,N_traj);

for i_traj = 1:N_traj
    
    x_est = x0; P = P0;
    x = mvnrnd(x_est', P0, 1)';
    
    X(:,1,i_traj) = x;
    z_bias_x = 2*(rand - 0.5)*bhw_x;
    z_bias_y = 2*(rand - 0.5)*bhw_y;
    
    for i = 2:N
        
        err = x_est - x_nom(:,i-1);
        u = u_nom - K*err;
        
        w = mvnrnd([0 0], Q, 1)';
        x = A*x + B*u + w;
        X(:,i,i_traj) = x;
        
        if mm_available
            if rand < p_bias_jump
                z_bias_x = 2*(rand - 0.5)*bhw_x;
                z_bias_y = 2*(rand - 0.5)*bhw_y;
            end
            z = C*x + [z_bias_x;z_bias_y] + mvnrnd([0 0], diag([z_cov_x z_cov_y]), 1)';
            Z(:,i,i_traj) = z;
            
            x_pred = A*x_est + B*u;
            P_pred = A*P*A' + Q;
            
            L = P_pred*C'/(C*P_pred*C' + R);
            x_est = x_pred + L*(z - C*x_pred);
            P = P_pred - L*C*P_pred;
        else            
            x_pred = A*x_est + B*u;
            P_pred = A*P*A' + Q;
            
            x_est = x_pred;
            P = P_pred;
        end

    end
    
end

%% RRBT calculation
G = zeros(2);
RRBT = cell(1,N);
P = P0;
RRBT{1} = P0;

for i = 2:N
        
    if mm_available
        P_pred = A*P*A' + Q;
        L = P_pred*C'/(C*P_pred*C' + R);
        P = P_pred - L*C*P_pred;
        G = (A-B*K)*G*(A-B*K)' + L*C*P_pred;
    else
        P = A*P*A' + Q;
        G = (A-B*K)*G*(A-B*K)';
    end
    
    RRBT{i} = P+G;
    
end



%% plot
m = 3;

RRBT_ms = nan(1,N);
for i = 1:N
    cov = RRBT{i};
    RRBT_ms(1,i) = m*sqrt(cov(2,2));
end

figure()
hold on; grid on;
for i_traj = 1:N_traj
    plot(X(1,:,i_traj),X(2,:,i_traj),'b');
end
plot(1:N,RRBT_ms,'r','Linewidth',2)
plot(1:N,-RRBT_ms,'r','Linewidth',2)
axis equal
xlabel('x-coordinate (m)');
ylabel('y-coordinate (m)');
