% this function search each S indices from T and gives matching T indicies

function T_ind=find_index_bh(T,S)

for i=1:length(S)
    try

T_ind(i)=find(T==S(i));
    catch
        if ~isempty(find(T==S(i)))
        error('Search vector is not unique')
        else
        T_ind(i)=0;    
        end
end
end
end