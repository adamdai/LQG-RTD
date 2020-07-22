function [] = generate_linear_dynamics()
% dynamics -- x and y position, heading, velocity, acceleration, are the
% states, and acceleration and yaw rate (both assumed constant over this
% time horizon) are the trajectory parameters.

   syms x y vx vy; % variables for FRS
   dx = vx;
   dy = vy;
   
   z = [x; y];
   dz = [dx; dy];
   
   syms tdummy udummy % dummy variables
   
   matlabFunction(dz, 'File', 'dyn_lin', 'vars', {tdummy z udummy});

end