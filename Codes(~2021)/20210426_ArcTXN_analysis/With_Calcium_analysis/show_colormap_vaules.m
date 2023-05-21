function show_colormap_vaules(values,min_max)
cmap=jet(100);
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
            'Color',[1 1 1],...
            'Renderer','painters','position',[100 100 200 200]);
        if min(values)<=min_max(1)
            min_max(1)=min(values)*0.9;
        end
        values(find(values>min_max(2)))=min_max(2);
for v=1:length(values)

plot(v,1,'marker','.','color',cmap(ceil((values(v)-min_max(1,1))/(min_max(1,2)-min_max(1,1))*100),:),'markersize',50)
hold all
end
axis equal tight off
end