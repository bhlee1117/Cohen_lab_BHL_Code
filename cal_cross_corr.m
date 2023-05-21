function CorrSpatial_vectors=cal_cross_corr(traces,coord,threshold)
%CorrSpatial_vectors=zeros(size(traces,1)*(size(traces,1)-1)/2,4);
g=1;
for i=1:size(traces,1)
    for j=i+1:size(traces,1)
    corrMat(i,j,:)=xcorr(traces(i,:),traces(j,:),300,'biased');
    [m(i,j) arg]=max(corrMat(i,j,:),[],3);
    dist_xy(i,j,1)=coord(i,1)-coord(j,1);
    dist_xy(i,j,2)=coord(i,2)-coord(j,2);
    if m(i,j)>threshold
    CorrSpatial_vectors(g,:)=[coord(i,:) -squeeze(dist_xy(i,j,1:2))'*sign(arg)];
    g=g+1;
    end
    
    end
end
end
