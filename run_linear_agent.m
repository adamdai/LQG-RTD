% Parameterize trajectory space and sample a few trajectories
% and track them

%A = linear_agent();
A = discrete_agent();
x0 = [0;0];
A.reset(x0)

t_f = 1;
vx_des = 1.0; 
vy_des = 2.0; 
[T,U,Z] = make_linear_trajectory(t_f,vx_des,vy_des);

t_total = t_f;
A.move(t_total,T,U,Z)

%% plotting
figure(1) ; clf ; axis equal ; hold on ; set(gca,'FontSize',15)

plot(Z(1,:),Z(2,:),'b--','LineWidth',1.5)
plot(A)
A.animate()






