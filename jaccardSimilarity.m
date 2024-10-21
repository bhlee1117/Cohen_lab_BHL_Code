function jaccardIndex = jaccardSimilarity(cluster1, cluster2)
    % Ensure both cluster assignments are of the same size
    assert(length(cluster1) == length(cluster2), 'Cluster assignments must have the same length');
    
    % Compute pairs of points that are clustered together in both clusterings
    n = length(cluster1);
    pairsSameCluster1 = zeros(n, n);
    pairsSameCluster2 = zeros(n, n);
    
    for i = 1:n
        for j = 1:n
            pairsSameCluster1(i, j) = (cluster1(i) == cluster1(j));
            pairsSameCluster2(i, j) = (cluster2(i) == cluster2(j));
        end
    end
    
    % Jaccard Index: intersection over union
    intersection = sum(sum(pairsSameCluster1 & pairsSameCluster2));
    union = sum(sum(pairsSameCluster1 | pairsSameCluster2));
    
    jaccardIndex = intersection / union;
end
