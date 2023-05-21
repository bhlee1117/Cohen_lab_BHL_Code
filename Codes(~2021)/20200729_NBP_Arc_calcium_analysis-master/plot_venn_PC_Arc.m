%function plot_venn_PC_Arc
xtick={'~PC & ~Arc','~PC & Arc','PC & ~Arc','PC & Arc'};


for m=1:size(dat,2)
    Arc_class=zeros(size(dat{m}.Arc,1),3);
    Arc_class(dat{m}.Arc(:,2:end)==5)=NaN;
    Arc_class(dat{m}.Arc(:,2:end)==1 | dat{m}.Arc(:,2:end)==3)=0;
    Arc_class(dat{m}.Arc(:,2:end)==2 | dat{m}.Arc(:,2:end)==4)=1;
    PC=cell2mat(dat{m}.isPC);
    for d1=1:3
    overlap{m}(d1,:)=[sum(Arc_class(:,d1)==0 & PC(:,d1)==0) ...
                      sum(Arc_class(:,d1)==1 & PC(:,d1)==0) ...
                      sum(Arc_class(:,d1)==0 & PC(:,d1)==1) ...
                      sum(Arc_class(:,d1)==1 & PC(:,d1)==1)]...
                    ./sum(~isnan(Arc_class(:,d1)) & ~isnan(PC(:,d1)));
    end
end
%%
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 300 300]);
Ov=cell2mat(overlap');
M=mean(Ov([1:12 16 17 18],:));
S=std(Ov([1:12 16 17 18],:));
cmap=distinguishable_colors(6);
[a Venn_data]=venn(M(2:end));
hold all
text(-1,-1,num2str(round(M(1,1)*100)/100,1),'HorizontalAlignment', 'center','Fontsize',11,'FontName','arial rounded mt bold')
for i=2:size(M,2)
    text(Venn_data.ZoneCentroid(i-1,1),Venn_data.ZoneCentroid(i-1,2),...
    num2str(round(M(1,i)*100)/100),'HorizontalAlignment', 'center','Fontsize',11,'FontName','arial rounded mt bold')
end
axis equal off tight
figure2 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 300 300]);
errorbar([1:size(M,2)],M,S,'LineWidth',2,'linestyle','none','color','k','Capsize',10)
hold all
for i=1:size(Ov,1)
    plot([1:size(M,2)],Ov(i,:),'color',cmap(ceil(i/3),:),'marker','.','linewidth',2)
    hold all
end
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:size(M,2)],'XTickLabel',...
 xtick,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
ylabel('Fraction of neurons','LineWidth',2,'FontSize',13,...
    'FontName','arial rounded mt bold')
xlim([0.5 size(M,2)+0.5])
ylim([0 1])
xtickangle(45)
%end