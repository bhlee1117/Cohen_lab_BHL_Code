function [identified_list identified_calcium identified_SpC identified_spike identified_calcium_raw]=identify_sp_comp(result,match_CS)
g=[];
for i=1:size(match_CS,2)
[a order]=sort(match_CS{1,i}(:,1),'ascend');
as_match_CS{1,i}=match_CS{1,i}(order,:);
g=[g;unique(a)];
end
g=unique(g);
for i=1:size(g,1) %Cell ID
    identified_list(i,1)=g(i,1);
    for j=1:size(match_CS,2) %Number of days
        if isempty(find(match_CS{1,j}(:,1)==g(i,1)))
        identified_list(i,j+1)=NaN;    
        identified_calcium{i,j}=[];
        else
        identified_list(i,j+1)=match_CS{1,j}(find(match_CS{1,j}(:,1)==g(i,1)),2);
        identified_calcium{i,j}=result{j}.C_df(match_CS{1,j}(find(match_CS{1,j}(:,1)==g(i,1)),2),:);
        identified_SpC{i,j}=result{j}.A_or(:,match_CS{1,j}(find(match_CS{1,j}(:,1)==g(i,1)),2));
         identified_spike{i,j}=   result{j}.S_or(match_CS{1,j}(find(match_CS{1,j}(:,1)==g(i,1)),2),:);
         identified_calcium_raw{i,j}=   result{j}.S_or(match_CS{1,j}(find(match_CS{1,j}(:,1)==g(i,1)),2),:);
        end
    end
end
end