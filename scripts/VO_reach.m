%% Standard VO

% read in images
I1 = imread('000000.png');
I2 = imread('000001.png');

% extract features
F1 = detectFASTFeatures(I1);
F2 = detectFASTFeatures(I2);

figure(1);
imshow(I1); hold on;
plot(F1.selectStrongest(50));

figure(2);
imshow(I2); hold on;
plot(F2.selectStrongest(50));

% perform 


%% VO reach

% imread img 1
% imread img 2

% add noise to imgs and convert to zonotopes 

% extract features