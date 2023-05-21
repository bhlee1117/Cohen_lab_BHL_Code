function coord=get_coord(c_ftprnt)

[x y]=meshgrid([1:size(c_ftprnt,2)],[1:size(c_ftprnt,1)]);
for i=1:size(c_ftprnt,3)
%coord(i,:)=[mean(find(max(c_ftprnt(:,:,i),[],1)>0)) mean(find(max(c_ftprnt(:,:,i),[],2)>0))];
coord(i,:)=[sum(sum(c_ftprnt(:,:,i).*x))./sum(sum(c_ftprnt(:,:,i))) sum(sum(c_ftprnt(:,:,i).*y))./sum(sum(c_ftprnt(:,:,i)))];
end
end