function MovM=ringmovMean(M,g)

MovM=repmat(M,1,3);
MovM=movmean(MovM,g,2,'omitnan');
% gF = fspecial('gaussian', 7, floor(g/2));
% MovM = imfilter(MovM, gF(4,:), 'same', 'replicate');
MovM=MovM(:,size(M,2)+1:2*size(M,2),:);

end