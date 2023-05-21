function [ident_list match]=match_cell(C,d_th)

if size(C,2)==1
ident_list=[1:size(C{1,1},3)]';
match=[];
else

for i=1:length(C)
coord{i}=get_coord(C{i});
P{i}=ones(size(coord{i},1),1);
end
ind=[1:length(C)];
pairs=nchoosek(ind, 2);

for i=1:size(pairs,1)
    n=[size(coord{pairs(i,1)},1) size(coord{pairs(i,2)},1)];
    D=squareform(pdist([coord{pairs(i,1)}; coord{pairs(i,2)}]));
    D=D(1:n,n+1:end); [a b]=find(D<d_th);
    r=[a b];
    match{i}=r;
    for j=1:2       
    t=zeros(n(j),1); t(r(:,j),1)=1;
    P{pairs(i,j)}=[P{pairs(i,j)}.*t];
    end
    for k=2:length(C)
    if sum(ismember(pairs(i,:),[1 k]))==2
    s(k)=i;
    end
    end
end

x=find(P{1})
ident_list=[x];
for k=2:length(C)
[a b]=(ismember(match{s(k)}(:,1),x));
ident_list=[ident_list match{s(k)}(b(find(a)),2)];
end


end
end