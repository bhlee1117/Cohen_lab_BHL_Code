function [match_list unmatched_list match_center]=coordinate_match(cell_list,center,max_dist,argz,range,im)
CL=cell_list(find(cell_list(:,3)>argz-range&cell_list(:,3)<argz+range),:);
CL_arg=find(cell_list(:,3)>argz-range&cell_list(:,3)<argz+range);
center=center*2;
for i=1:size(CL,1)
    for j=1:size(center,1)
        D(i,j)=distance_BH(CL(i,:),[center(j,[1 2]) argz]);
    end
end

ascend_dist=sort(D(:),'ascend');
match_list=nan(1,2);
g=1;
for i=1:size(find(ascend_dist<max_dist),1)
    [r c]=find(D==ascend_dist(i,1));
    if isempty(find(match_list(:,1)==CL_arg(r,1))) && isempty(find(match_list(:,2)==c))
        match_list(g,1:2)=[CL_arg(r,1) c];
        match_center{g,1}=CL(r,:);
        match_center{g,2}=[center(c,[1 2]) argz];
        g=g+1;
    end
end
unmatched_list{1}=setdiff(1:size(cell_list,1),match_list(:,1));
unmatched_list{2}=setdiff(1:size(center,1),match_list(:,2));
if ~isnan(match_list)
draw_match(cell_list,center(:,[2 1]),match_list,unmatched_list,argz,im)
end
end
function draw_match(cell_list,center,match_list,unmatched_list,argz,im)
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 700 900]);
if ~isempty(im)
imagesc(im)
axis equal tight off
hold all
end

plot3(cell_list(match_list(:,1),1),cell_list(match_list(:,1),2),cell_list(match_list(:,1),3),'marker','.','color',[0.1 0.45 1],'markersize',20,'linestyle','none')
hold all
axis equal
plot3(cell_list(unmatched_list{1},1),cell_list(unmatched_list{1},2),cell_list(unmatched_list{1},3),'marker','.','color',[0.1 0.1 .3],'markersize',20,'linestyle','none')
plot3(center(match_list(:,2),2),center(match_list(:,2),1),repmat(argz,size(match_list,1),1),'marker','o','color',[1 0.45 1],'linewidth',2,'linestyle','none');
plot3(center(unmatched_list{2},2),center(unmatched_list{2},1),repmat(argz,size(unmatched_list{2},2),1),'marker','o','color',[0.3 0.15 0.1],'linewidth',2,'linestyle','none');
for i=1:size(match_list,1)
    line([cell_list(match_list(i,1),1) center(match_list(i,2),2)],[cell_list(match_list(i,1),2) center(match_list(i,2),1)],[cell_list(match_list(i,1),3) argz],...
          'color',[0.6 0.2 0.6],'linewidth',2);
      hold all
end
grid on
end