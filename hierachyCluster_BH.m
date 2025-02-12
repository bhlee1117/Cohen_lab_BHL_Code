function [ExcPatt_reduce num_clusters clusterBinMat]= hierachyCluster_BH(ExcPatt,cuttoff_val,figureSwitch)

if nargin<3
    figureSwitch=0;
end
distances = pdist(ExcPatt', 'euclidean');  % Euclidean distance
Z = linkage(distances, 'ward');
%num_clusters = 15;
%cuttoff_val=24;
cluster_indices = cluster(Z, 'MAXCLUST', cuttoff_val, 'depth', 2);
sortedHeights = sort(Z(:, 3), 'descend'); % Linkage heights
cutoff = sortedHeights(cuttoff_val - 1);

num_clusters=max(cluster_indices);

ExcPatt_reduce = zeros(size(ExcPatt, 1), num_clusters); cluster_weight=zeros(1,num_clusters);
for i = 1:num_clusters
    ExcPatt_reduce(:, i) = mean(ExcPatt(:, find(cluster_indices == i)),2,'omitnan');
    %Vpre_reduce(:, i) = Vpre(:, find(cluster_indices == i,1));
    cluster_weight(i)=length(find(cluster_indices == i));
end
[~,sort_cluster]=sort(cluster_weight,'descend');
ExcPatt_reduce=ExcPatt_reduce(:,sort_cluster);
[VExcpatt, ~, ExcPattTrace]=get_eigvector(ExcPatt');
clusterBinMat=get_indMat(find_index_bh(sort_cluster,cluster_indices)');

% Alternatively, you can plot a dendrogram to visualize the clustering
%imagesc(rescale2(Vs_reduce(:,sort_cluster),1))
if figureSwitch
clf; tiledlayout(1,3);
nexttile([1 1]); dendrogram(Z,size(Z,1),'ColorThreshold',cutoff);
nexttile([1 1]); silhouette(ExcPatt',cluster_indices,'Euclidean');
nexttile([1 1]); 
for c=1:num_clusters
    scatter3(ExcPattTrace(find(clusterBinMat(:,c)),1),ExcPattTrace(find(clusterBinMat(:,c)),2),ExcPattTrace(find(clusterBinMat(:,c)),3),'filled'); hold all
end
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
end
end