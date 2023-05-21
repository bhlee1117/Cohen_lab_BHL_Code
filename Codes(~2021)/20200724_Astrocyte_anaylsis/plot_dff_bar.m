function p_value=plot_dff_bar(DFF,stim_time,timescale,xticks,color)
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 400 300]);
ylim=5;
grid=[0.8 1.2 1.8 2.2];
for i=1:size(DFF,1)
    for j=1:size(DFF,2)
datum{1}(j,i)=mean(mean(DFF{i,j}(:,1:round(stim_time/timescale)))); % Basal
datum{2}(j,i)=mean(mean(DFF{i,j}(:,round(stim_time/timescale):end))); % After stimulation
    end
end
M=[mean(datum{1}) mean(datum{2})];
S=[std(datum{1},0,1)/sqrt(size(datum{1},1)) std(datum{2},0,1)/sqrt(size(datum{2},1))];
errorbar(grid,M,S,'LineWidth',2,'linestyle','none','color','k','Capsize',10)
hold all
b=bar(grid,M,'Barwidth',0.7,'LineWidth',2);
grid
set(gca,'FontSize',8,'LineWidth',1,'XTick',grid,'XTickLabel',...
 xticks,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
b.FaceColor='flat';
for j=1:4
b.CData(j,:)= color(j,:);
end
for i=1:2
[a p_value(i,1)]=ttest2(datum{1}(:,1),datum{1}(:,2));
end
ref=1;
sh=0.2;
[r c]=find(p_value(:,ref)<0.01);
[rr cc]=find(p_value(:,ref)<0.05 & p_value(:,ref)>0.01);
for i=r'
    line([i-sh i+sh],[y_lim*(0.87+0.1*(ref-1)) y_lim*(0.87+0.1*(ref-1))],'color','k','LineWidth',1.5);
    line([i-sh i-sh],[100*(M(1,2*i-1)+S(1,2*i-1))+3 y_lim*(0.87+0.1*(ref-1))],'color','k','LineWidth',1.5);
    line([i+sh i+sh],[100*(M(1,2*i)+S(1,2*i))+3 y_lim*(0.87+0.1*(ref-1))],'color','k','LineWidth',1.5);
   text(i,y_lim*(0.89+0.1*(ref-1)),'**','HorizontalAlignment', 'center'...
         ,'LineWidth',2,'FontSize',13,'FontName','arial rounded mt bold')
end
for i=rr'
    line([i-sh i+sh],[y_lim*(0.87+0.1*(ref-1)) y_lim*(0.87+0.1*(ref-1))],'color','k','LineWidth',1.5);
    line([i-sh i-sh],[100*(M(1,2*i-1)+S(1,2*i-1))+3 y_lim*(0.87+0.1*(ref-1))],'color','k','LineWidth',1.5);
    line([i+sh i+sh],[100*(M(1,2*i)+S(1,2*i))+3 y_lim*(0.87+0.1*(ref-1))],'color','k','LineWidth',1.5);
    text(i,y_lim*(0.89+0.1*(ref-1)),'*','HorizontalAlignment', 'center'...
         ,'LineWidth',2,'FontSize',13,'FontName','arial rounded mt bold')
end


ylabel('Mean \DeltaF/F','FontSize',13,'FontName','arial rounded mt bold','HorizontalAlignment', 'center')
xlabel('     Basal                 Encounter','FontSize',13,'FontName','arial rounded mt bold','HorizontalAlignment', 'center')

end