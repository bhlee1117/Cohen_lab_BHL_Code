function match_cell_spatial_comp(cell_list,result,z_cal,options,min_dist) 
cal_centroid=2*cell2mat(cellfun(@mean_BH,result.Coor','UniformOutput',false))';
for i=1:size(cell_list,1) %Arc_TXN cell
    for j=1:size(cal_centroid,1) %Spatial component
        distance(i,j)=distance_BH(cell_list(i,:),[cal_centroid(j,:) z_cal]);
    end
end
[closest matched]=min(distance,[],2);
match_list=find(closest<=min_dist);
unmatched_list=setdiff(1:size(cal_centroid,1),matched(match_list,1));
spatial_comp_plot(result.A_or,options,[],2) % color, scale
hold all
plot(cell_list(match_list,1),cell_list(match_list,2),'marker','.','color',[0.1 0.45 1],'markersize',20,'linestyle','none')
plot(cell_list(find(closest>min_dist),1),cell_list(find(closest>min_dist),2),'marker','.','color',[0.4 0.4 0.4],'markersize',20,'linestyle','none')
plot(cal_centroid(matched(match_list,1),1),cal_centroid(matched(match_list,1),2),'marker','o','color',[1 0.45 1],'linewidth',2,'linestyle','none');
plot(cal_centroid(unmatched_list,1),cal_centroid(unmatched_list,2),'marker','o','color',[0.5 0.35 0.5],'linewidth',2,'linestyle','none');
for i=1:size(match_list,1)
    line([cell_list(match_list(i,1),1) cal_centroid(matched(match_list(i,1),1),1)],[cell_list(match_list(i,1),2) cal_centroid(matched(match_list(i,1),1),2)],...
          'color',[0.6 0.2 0.6],'linewidth',2);
      hold all
end
end