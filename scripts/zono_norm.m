%% description
% test find the interval imposed by taking the norm over all the points in 
% a zonotope

% c = [3;3];
% G = [1 0.5;
%      0.2 2];
% Z = zonotope(c,G);
% 
% plot(Z)

function I = zono_norm(Z)

% sample points
N = 1000;
X = sampleBox(Z,N);
min_norm = inf;
max_norm = 0;
for i = 1:N
    x = X(:,i);
    d = norm(x);
    if d < min_norm
        min_norm = d;
    end
    if d > max_norm
        max_norm = d;
    end
end

I = interval(min_norm, max_norm);

end

