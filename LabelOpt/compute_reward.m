function r = compute_reward(points)
% Computes the reward (sum of pairwise Euclidean distances)
try
    distances = pdist(points, 'euclidean'); % compute pairwise distances vectorized
    r = sum(distances); % sum all distances
catch
    warning('Vectorized computation failed (likely due to memory constraints). Using loops instead.');
    r = 0;
    n = size(points, 1);
    for t1 = 1:n-1
        for t2 = t1+1:n
            r = r + norm(points(t1,:) - points(t2,:));
        end
    end
end
end
