% trajectory discretization and length
t_f = 5; dt = 0.01; N = t_f/dt; 

% system
A = eye(2); B = dt*eye(2); C = eye(2);
K = dlqr(A,B,eye(2),eye(2)); % feedback law u = -Kx
Q = 1e-5*eye(2); % process noise covariance
R = 1e-5*eye(2); % measurement noise covariance

% form sys struct
sys.A = A; sys.B = B; sys.C = C;
sys.K = K; sys.Q = Q; sys.R = R;

x0 = [0;0]; % inital state
u = [0.5; 0.5]; % control input

N_traj = 500;

% test for different initial covariances

% %-------------------------------------------------------------------------%
% P0 = 1e-3*eye(2); 
% X = simulate_LQG_trajectory(sys,N,N_traj,u,x0,P0);
% 
% figure(1);
% hold on; grid on;
% for i_traj = 1:N_traj
%     plot(X(1,:,i_traj),X(2,:,i_traj),'r');
% end
% axis equal
% xlabel('x-coordinate (m)');
% ylabel('y-coordinate (m)');
% 
% %-------------------------------------------------------------------------%
% P0 = 1e-2*eye(2); 
% X = simulate_LQG_trajectory(sys,N,N_traj,u,x0,P0);
% 
% figure(2);
% hold on; grid on;
% for i_traj = 1:N_traj
%     plot(X(1,:,i_traj),X(2,:,i_traj),'b');
% end
% axis equal
% xlabel('x-coordinate (m)');
% ylabel('y-coordinate (m)');

%-------------------------------------------------------------------------%
P0 = 1e-3*diag([1 10]); 
X = simulate_LQG_trajectory(sys,N,N_traj,u,x0,P0);

figure(3);
hold on; grid on;
for i_traj = 1:N_traj
    plot(X(1,:,i_traj),X(2,:,i_traj),'b');
end
axis equal
xlabel('x-coordinate (m)');
ylabel('y-coordinate (m)');