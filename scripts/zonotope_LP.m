%% parameters


%% Standard LP
% max_{x}     c^T x
% subject to  Ax <= b
%             x >= 0

c = [0.6, 0.35];
A = [5 7;
     4 2;
     2 1];
b = [8; 15; 3];


% simplex method (?)