function [ref_result X]=find_baseline(result)

for i=1:size(result.C_df,1)
plot([1:1:size(result.C_df(i,:),2)],result.C_df(i,:))
[x{i} y]=ginput(2);
end
g=1;
for i=1:size(x,2)
if isempty(x{i})
disp([num2str(i) ' th cell was removed'])
else
ref_result.C_df(g,:)=result.C_df(i,:);
ref_result.A_or(:,g)=result.A_or(:,i);
ref_result.C_or(g,:)=result.C_or(i,:);
ref_result.S_or(g,:)=result.S_or(i,:);
ref_result.Coor{g,1}=result.Coor{i,1};

X(g,:)=x{i};
g=g+1;
end
end
ref_result.options=result.options;
end