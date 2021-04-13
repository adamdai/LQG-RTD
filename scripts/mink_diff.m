%% description
% test cases of minkowski difference

%% CORA example

Z1 = zonotope([1 2 2; 0 0 2]);
Z2 = zonotope([0 0.5 0.5 0.3; 1 0 0.5 0.2]);
Z3 = Z1 - Z2;
Z4 = Z2 + Z3;
plot(Z1,[1 2], 'b');
hold on
plot(Z2,[1 2], 'r');
plot(Z3,[1 2], 'g');
plot(Z4,[1 2], 'k');

%% generic example

Z1 = zonotope([0 1 0; 
               0 0 1]);
Z2 = zonotope([1.5 1   0.5; 
               1.5 0.5   1]);
Z3 = Z1 + Z2;

Z4 = Z3 - Z1; 
Z5 = Z3 - Z2; 

figure(1); hold on; axis equal
plot(Z1,[1,2],'FaceColor','r','FaceAlpha',0.1,'Filled',true);
plot(Z2,[1,2],'FaceColor','g','FaceAlpha',0.1,'Filled',true);
plot(Z3,[1,2],'FaceColor','b','FaceAlpha',0.1,'Filled',true);

figure(2); hold on; axis equal
plot(Z4,[1,2],'FaceColor','r','FaceAlpha',0.1,'Filled',true);
plot(Z5,[1,2],'FaceColor','g','FaceAlpha',0.1,'Filled',true);

Z6 = Z5 + Z2;
plot(Z6,[1,2],'FaceColor','b','FaceAlpha',0.1,'Filled',true);

figure(3); hold on; axis equal
Z7 = Z2 + Z1;
plot(Z3,[1,2],'FaceColor','r','FaceAlpha',0.1,'Filled',true);
plot(Z7,[1,2],'FaceColor','g','FaceAlpha',0.1,'Filled',true);

%% error case

c_g = [0; 0];
G_g = [1; 0];
g = zonotope(c_g,G_g);

c_y = [0; 0];
G_y = eye(2);
y = zonotope(c_y,G_y);

g-y

hold on; grid on; axis equal
plot(g,[1,2],'FaceColor','r','FaceAlpha',0.1,'Filled',true);
plot(y,[1,2],'FaceColor','b','FaceAlpha',0.1,'Filled',true);
