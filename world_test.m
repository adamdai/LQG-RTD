%% description
% This script runs a simulation with an LQG agent in the simulator
% framework, using RTD to plan online.
%
% Author: Adam Dai
% Created: 30 July 2020
%
%% user parameters
% world
obstacle_size_bounds = [0.5, 0.5] ; % side length [min, max]
N_obstacles = 7 ;
bounds = [-4,4,-2,2] ;
goal_radius = 0.5 ;

% planner
buffer = 0.05 ; % m
t_plan = 0.5 ; % if t_plan = t_move, then real time planning is enforced
t_move = 0.5 ;

% simulation
verbose_level = 10 ;

%% automated from here
A = LQG_agent ;

W = static_box_world_noheading('bounds',bounds,'N_obstacles',N_obstacles,'buffer',0.25,...
                     'verbose',verbose_level,'goal_radius',goal_radius,...
                     'obstacle_size_bounds',obstacle_size_bounds) ;

figure(1);
cla ; hold on ; axis equal ;
axis(W.bounds)
plot(W)