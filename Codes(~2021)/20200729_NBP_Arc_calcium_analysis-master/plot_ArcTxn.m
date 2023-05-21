function [on_frac2 Fraction]=plot_ArcTxn(Arcclass,post_ref,cmap,y_lim,xtick) %sw for line
figure2 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 350 350]);
G=[0.8 1.2 1.8 2.2 2.8 3.2];
for m=1:size(Arcclass,2)
    for d1=1:3
        for k=1:4
        Fraction{1,m}(k,d1)=sum(Arcclass{m}(:,d1)==k)/size(Arcclass{m},1);
        end
        Fraction{2,m}=[sum(Fraction{1,m}([1 3],:));sum(Fraction{1,m}([2 4],:))];
    end
end

if post_ref    
else
mat_fraction=cell2mat(Fraction(1,:));
on_frac2=[mat_fraction(3,:);sum(mat_fraction([2 4],:),1)];
S=std(reshape(on_frac2,2,3,size(Arcclass,2)),0,3)/sqrt(6);
on_frac=mean(reshape(on_frac2,2,3,size(Arcclass,2)),3);

errorbar(G,on_frac(:)*100,S(:)*100,'LineWidth',2,'linestyle','none','color','k','Capsize',10)
hold all
b=bar(G,on_frac(:)*100,'Barwidth',0.7,'LineWidth',2);
b.FaceColor='flat';
for j=1:size(G,2)
b.CData(j,:)= cmap(j,:);
end    
end
set(gca,'FontSize',8,'LineWidth',1,'XTick',G,'XTickLabel',...
 xtick,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
xlim([0.5 3+0.5])
ylim([0 y_lim])
xtickangle(45)
ylabel('Fraction of Arc^+ neurons(%)','FontSize',13,'FontName','arial rounded mt bold')
end