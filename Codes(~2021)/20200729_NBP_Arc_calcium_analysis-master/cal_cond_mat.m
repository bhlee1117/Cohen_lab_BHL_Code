function [cond_matrix chan_matrix omitna_Arcclass]=cal_cond_mat(Full_result,on_off_cond)
 comb=[1 1;1 2;2 1;2 2];

    for m=1:size(Full_result,2) % Mouse
omit_list=find(sum(Full_result{m}.Arc_class(:,4:end)==3,2)==0);
omitna_class=[Full_result{m}.Arc_class(omit_list,4:end)];
day_end=(size(Full_result{m}.Arc_class,2)-3)/2;
omitna_Arcclass{m}=zeros(size(omitna_class,1),day_end);
for i=1:day_end % Day of interest
    for k=1:size(comb,1) % Classification in Arc TXN
        l=find(omitna_class(:,2*i-1)==comb(k,1) &...
               omitna_class(:,2*i)==comb(k,2));    % find the list (identified list standard) of each combination
        omitna_Arcclass{m}(l,i)=k;
    end
end
  omitna_Arcclass{m}(find(sum(omitna_Arcclass{m}==0,2)~=0),:)=[]; % Na
  clear l l2
 for j=1:day_end %overlap case
     for jj=1:day_end
             for sw=1:size(on_off_cond,2)
             l(:,sw)=omitna_Arcclass{m}(:,j) ==on_off_cond(1,sw);
            l2(:,sw)=omitna_Arcclass{m}(:,jj)==on_off_cond(1,sw); end
         l=sum(l,2)>0; l2=sum(l2,2)>0;    
         cond_matrix{1}{1,m}(j,jj)=sum(l&l2)/size(omitna_Arcclass{m},1); % overlap
         cond_matrix{2}{1,m}(j,jj)=sum(l&l2)/sum(l|l2); % similarity, 'Jaccard'
         cond_matrix{3}{1,m}(j,jj)=sum(l&l2)/sum(l); % cond_prob
         
         chan_matrix{1}{1,m}(j,jj)=sum(l)*sum(l2)/(size(omitna_Arcclass{m},1)^2);
         chan_matrix{2}{1,m}(j,jj)=sum(l2)/size(omitna_Arcclass{m},1); % cond chance
         chan_matrix{3}{1,m}(j,jj)=sum(~l & l2)/sum(~l); % cond negative
%          cond_matrix{5}{i,m}(j,jj)=acos(overlap_matrix(j,jj)/sqrt(sum(omitna_class(:,j+1)==2)*sum(omitna_class(:,jj+1)==2)))*180/pi; %Angle
%          cond_matrix{6}{i,m}(j,jj)=sum((omitna_class(:,j+1)==2) ~= (omitna_class(:,jj+1)==2))/size(omitna_class,1); %distance
         if j==jj && isnan(cond_matrix{3}{1,m}(j,jj))
             cond_matrix{3}{1,m}(j,jj)=1;
             chan_matrix{3}{1,m}(j,jj)=1;
         end
     end
 end
    end
end