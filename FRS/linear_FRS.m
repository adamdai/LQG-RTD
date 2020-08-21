% Given a nominal input range K1 x K2, compute the FRS for a 2D linear 
% trajectory model
% 
% Inputs:
%   K1 - +/- interval for k1 trajectory parameter (x-speed)
%   K2 - +/- interval for k2 trajectory parameter (y-speed)
%
% Outputs:
%   frs - cell array of zonotopes representing forward reachable sets in time
%   

function frs = linear_FRS(K1,K2,t_f,dt)

    % model initial condition set:
    % X0 = [x_0; % initial x position
    %       y_0; % initial y position
    %       K1;  % desired x speed (K1 parameter set, centered at k1c with generator k1g)
    %       K2]; % desired y speed (K2 parameter set, centered at k2c with generator k2g)

    % non-sliceable initial conditions
    x_0 = 0;
    y_0 = 0;

    % sliceable parameters... desired x speed and y speed
    k1c = 0; % center of x speed interval
    k1g = K1; % k1c +- k1g
    k2c = 0; % center of y speed interval
    k2g = K2; % k2c +- k2g

    % set FRS options
    options.tStart = 0;
    options.tFinal = t_f;
    options.x0 = [x_0; y_0; k1c; k2c]; % center of initial set
    options.R0 = zonotope([options.x0, [0;0;k1g;0], [0;0;0;k2g]]); % generators for parameter dimensions
    options.timeStep = dt; %time step size for reachable set computation
    options.taylorTerms = 5; %number of taylor terms for reachable sets
    options.zonotopeOrder = 20; %zonotope order... increase this for more complicated systems.
    options.maxError = 1*ones(4, 1); % this controls splitting, set it high to avoid splitting
    options.verbose = 1;

    options.uTrans = 0; % center of input set
    options.U = zonotope([0, 0]);

    options.originContained = 1;
    options.advancedLinErrorComp = 0;
    options.tensorOrder = 1;
    options.reductionInterval = inf;
    options.reductionTechnique = 'girard';

    % specify continuous dynamics----------------------------------------------
    % transition matrix for state augmented with constant input (k1,k2)
    A = [0 0 1 0; 
         0 0 0 1; 
         0 0 0 0; 
         0 0 0 0]; 
    B = zeros(4,1);
    sys = linearSys('sys',A,B);

    % compute reachable set----------------------------------------------------
    tic
    frs = reach(sys, options);
    tComp = toc;
    disp(['computation time of reachable set: ', num2str(tComp)]);
    
    % remove zero generators
    for i = 1: length(frs)
        frs{i} = deleteZeros(frs{i});
    end
    
end

