function data=venn_data(mat)
% mat is n x 3 matrix
data(1,1)=sum(mat(:,1)==2 & mat(:,2)==1 & mat(:,3)==1);
data(1,2)=sum(mat(:,1)==1 & mat(:,2)==2 & mat(:,3)==1);
data(1,3)=sum(mat(:,1)==1 & mat(:,2)==1 & mat(:,3)==2);
data(1,4:6)=[sum(mat(:,1)==2 & mat(:,2)==2 & mat(:,3)==1) sum(mat(:,1)==2 & mat(:,3)==2 & mat(:,2)==1) sum(mat(:,2)==2 & mat(:,3)==2 & mat(:,1)==1) ];
data(1,7)=sum(mat(:,1)==2 & mat(:,2)==2 & mat(:,3)==2);
end