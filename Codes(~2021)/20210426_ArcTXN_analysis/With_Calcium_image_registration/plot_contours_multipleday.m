function plot_contours_multipleday(identified_SpC,options,cell_list,identified_list,argz,doi)
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 1200 900]);
[a aa]=cellfun(@size,identified_SpC(:,doi));
%overlap_list=find(min(aa,[],2)>0)'; %Detected in All days
%cmap=distinguishable_colors(size(doi,2));
cmap=[0 0 1;1 0.6 0;1 0 0];
%cmap=[1 0 0];
for i=1:size(identified_SpC,1) %Cell
for j=doi %Day
    if ~isempty(identified_SpC{i,j})
A_temp = full(reshape(identified_SpC{i,j},options.d1,options.d2));
A_temp = medfilt2(A_temp,[3,3]);
A_temp = A_temp(:);
 [temp,ind] = sort(A_temp(:).^2,'ascend'); 
 temp =  cumsum(temp);
 ff = find(temp > (1-options.nrgthr)*temp(end),1,'first');
 p=full(reshape(identified_SpC{i,j},options.d1,options.d2));
 p(p > (1-options.nrgthr)*temp(end))=argz(j,1);
%  [s1 s2]=size(reshape(A_temp,result{j}.options.d1,result{j}.options.d2))
%   [s1 s2]=size([0,0]+A_temp(ind(ff)))
%  contour3(reshape(A_temp,result{j}.options.d1,result{j}.options.d2),[0,0]+A_temp(ind(ff)),zeros(433,300)+argz(j,1),'LineColor',cmap(j,:), 'linewidth', 2);
contour3(p,[argz(j,1) argz(j,1)],'LineColor',cmap(j,:), 'linewidth', 2);
 hold all
    end
end
%p3=plot3(cell_list(identified_list(i,1),1),cell_list(identified_list(i,1),2),cell_list(identified_list(i,1),3),'marker','.','markersize',50,'color',[0.3 0.3 0.3 0.1]);
set(gca,'FontSize',8,'LineWidth',1,'ZTick',[0:48:90],'ZTickLabel',...
 [0:48:90]*0.25,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2,...
'XTick',[0:60:700],'XTickLabel',[0:60:700]*0.25,'YTick',[0:60:1000],'YTickLabel',[0:60:1000]*0.25);
%  text(cell_list(identified_list(i,1),1),cell_list(identified_list(i,1),2),cell_list(identified_list(i,1),3)+10,num2str(identified_list(i,1)),'color','k'...
%  ,'FontName','arial rounded mt bold','Fontsize',14,'HorizontalAlignment', 'center')
end
zlim([0 90])
axis equal tight
end