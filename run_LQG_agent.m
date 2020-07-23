% Initialize an LQG agent (discrete or continuous) and simulate it tracking
% a nominal linear trajectory

rng(252)

A = LQG_agent();
%A = discrete_LQG_agent();
x0 = [0;0]; x_est0 = [0;0]; P0 = 5*eye(2);
A.reset([x0; x_est0; vec(P0)])

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






