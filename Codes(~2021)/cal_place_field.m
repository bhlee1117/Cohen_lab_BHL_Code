function [C_df cell_arranged]=cal_place_field(C_df,position,bin,min_track,range,method,arrange)

% g=1;
% for i=1:size(C_df,1)
% if isempty(find(C_df(i,:)>5))
% C_df_omit(g,:)=C_df(i,:);
% g=g+1;
% end
% end
if isempty(range)
    range=[1:size(C_df,1)];
end
switch method
    case 'Cdf'
        C_df=full(C_df([range],:));
    case 'NormCdf'
        %C_df=full((C_df-min(C_df,[],2))./(max(C_df,[],2)-min(C_df,[],2)));
        C_df=C_df([range],:);
end
place_cdf=[position(1:11000,1) C_df(:,1:11000)'];
bin_place=round((place_cdf(:,1)-min_track)/bin);
clear place_fr place_fr_gauss
for i=1:max(bin_place) %place
    for j=1:size(place_cdf,2)-1 %cell
    place_fr(i,j)=full(mean(place_cdf(find(bin_place==i),j+1)));
    end
end

for i=1:size(place_cdf,2)-1 %cell
    place_fr_gauss(i,:)=imgaussfilt(place_fr(:,i),1)';
end
place_fr_gauss=(place_fr_gauss-min(place_fr_gauss,[],2))./(max(place_fr_gauss,[],2)-min(place_fr_gauss,[],2));
[m max_FR_place]=max(place_fr_gauss,[],2);
[B cell_arranged]=sort(max_FR_place,'ascend');
 if arrange
imagesc(place_fr_gauss(cell_arranged,:),[0.2 1])
colormap('jet')
axis tight off
 else
imagesc(place_fr_gauss,[0.2 1])
colormap('jet')
axis tight off
 end
end
