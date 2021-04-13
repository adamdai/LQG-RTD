%% description
% This script tests a simple reachability example using constrained zonotopes.
% An initial constrained zonotope is created and propagated through a linear 
% dynamical system.
%
% Authors: Adam Dai
% Created: 5 Apr 2021
% Updated: 
%
%% user parameters

% initial constrained zonotope
% zonotope
Z = [0 1 0 2;
     0 0 1 -1]; 
% constraints
Z_A = [-2 1 -1];
Z_b = 2;
 
% system
A = [1 0.1;
     0.1 1];
B = [0.5 0;
     0 0.5];

% nominal trajectory
N = 10;
U = ones(2,N);

% norm to consider
p = 2 ;

x_0 = conZonotope(Z,Z_A,Z_b);

%% automated from here
% reachability
R = {};

x = x_0;
for i = 1:N
    u_nom = U(:,i);
    u = u_nom;
    
    % here, the * and + operators are overloaded with CORA linear mapping
    % and minkowski sum operators
    x = A * x + B * u; 
    
    R{i} = x;
end

%% plotting
hold on
for i = 1:N
   plot(R{i}); 
end