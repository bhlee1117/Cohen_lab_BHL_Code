function output=min_rank(mat,rank)
% mat is column vector
if size(mat,1)==1
    mat=mat';
end
tmp=mat;
for i=1:rank
    [m ind]=min(tmp);
    output(i,1)=m;
    tmp(ind,1)=max(mat);
end
end