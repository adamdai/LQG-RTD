% Parameterize trajectory space and sample a few trajectories
% and track them

A = turtlebot_agent();
v0 = 0.5;
A.reset([0;0;0;0.5])

% time horizon
t_f = 1; 
t_total = 1;

% trajectory space
K1 = -2:0.1:2; % yaw rate (rad/s)
K2 = 0:0.1:1; % speed (m/s)

figure(1) ; clf ; axis equal ; hold on ; set(gca,'FontSize',15)

for k1 = randsample(length(K1),2)'
    for k2 = randsample(length(K1),2)'
        [T,U,Z] = make_turtlebot_desired_trajectory(t_f,k1,k2);
        plot(Z(1,:),Z(2,:),'b--','LineWidth',1.5)
        
        A.move(t_total,T,U,Z)
        plot(A)
        A.reset([0;0;0;0.5])
    end
end

%% plotting
% figure(1) ; clf ; axis equal ; hold on ; set(gca,'FontSize',15)
% 
% plot(Z(1,:),Z(2,:),'b--','LineWidth',1.5)
% plot(A)
% A.animate()






