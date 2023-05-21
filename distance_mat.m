function dist=distance_mat(X,Y)
for i=1:size(X,1)
    for j=1:size(Y,1)
        dist(i,j)=sqrt(sum((X(i,:)-Y(j,:)).^2));
    end
end
end