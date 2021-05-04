%% parameters


%% Standard LP
% max_{x}     c^T x
% subject to  Ax <= b
%             x >= 0

A = [3,1,1,0;
     1,2,0,1];
b = [8,9];
c = [-4,-5,0,0];
format short;

% simplex method (?)
nma_simplex(A,b,c,false)