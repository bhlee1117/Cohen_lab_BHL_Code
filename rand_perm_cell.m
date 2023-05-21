function [perm_list perm]=rand_perm_cell(centers,N,dist_threshold)
i=1;
clist=centers;
perm_list=[];
perm=[];
while i<N
    if i>1
        n=randperm(size(clist,1),1);
        g=1;
        while sum(sqrt(sum((perm_list-clist(n,:)).^2,2))<dist_threshold)>0
            n=randperm(size(clist,1),1);
            
            g=g+1;
            if g>size(clist,1)
                i
                break;
            end
        end
        perm(i)=n;
    perm_list=[perm_list; clist(n,:)];
    clist(n,:)=[];
    if isempty(clist)
        break;
    end
    else
        size(clist,1)
    n=randperm(size(clist,1),1);
    perm_list=[clist(n,:)];
    perm(i)=n;
    end
    i=i+1;
end
end

