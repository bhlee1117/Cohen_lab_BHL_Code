function plot_cond_bar
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 250 250]);
int_g_col_day=[1 5 2;1 5 3];
cmap=distinguishable_colors(size(int_g_col_day,1));

switch int_g_col_day(1,2)
    case 1
        ylab='Sum of Peaks (\sigma)';
    case 2
        ylab='Sum of Peaks (\DeltaF/F)';
    case 3
        ylab='Mean of Peaks (\sigma)';
    case 4
        ylab='Mean of Peaks (\DeltaF/F)';
    case 5
        ylab='Ca^2^+ event rate (Hz)';
end
for i=1:size(int_g_col_day,1)
    dat{i}=result{int_g_col_day(i,1)}(:,5*(int_g_col_day(i,3)-1)+int_g_col_day(i,2));
    M(i,1)=mean(dat{i},'omitnan');
    S(i,1)=std(dat{i},0,1,'omitnan')/sqrt(size(dat{i},1));
end
G=[1:1:size(int_g_col_day,1)];
errorbar(G,M,S,...
        'LineWidth',2,'linestyle','none','color','k','Capsize',10)
    hold all
b=bar(G,M,'Barwidth',0.7,'LineWidth',2);
b.FaceColor='flat';
for i=1:size(int_g_col_day,1)
X=randn(size(dat{i},1),1)*0.01+G(i);
b.CData(i,:)= cmap(i,:);
hold all
plot(X,dat{i},'marker','.','markersize',10,'linestyle','none','color',cmap(i,:)/2)
end

set(gca,'FontSize',8,'LineWidth',1,'XTick',G,'XTickLabel',...
 xtick,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
ylabel(ylab,'LineWidth',2,'FontSize',13,...
    'FontName','arial rounded mt bold')
xlim([0.5 size(int_g_col_day,1)+ 0.5])
xtickangle(45)

