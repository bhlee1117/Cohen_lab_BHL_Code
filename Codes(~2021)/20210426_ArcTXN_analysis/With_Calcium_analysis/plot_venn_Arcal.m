function Venn_data=plot_venn_Arcal(Arcclass,group,on_cond)
xtick={'HC','CFC','Ret 1','Ret 2','Ret 3','Ret 4','Ret 5','Ret 6','Remote'};
% figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
%     'Color',[1 1 1],...
%     'Renderer','painters','position',[100 100 250 250]);

c=[2 1 1;1 2 1;1 1 2;2 2 1;2 1 2;1 2 2; 2 2 2];
nc=[1 2 1 1 1 1;1 1 1 2 1 1;1 1 1 1 1 2;1 2 1 2 1 1;1 2 1 1 1 2;1 1 1 2 1 2;1 2 1 2 1 2];
% Make NA cell omited data
for m=1:size(Arcclass,2)
for d1=1:3
    clear l
    for sw=1:size(on_cond,2)
        l(:,sw)=Arcclass{m}(:,d1)==on_cond(1,sw);
    end
    l2{m}(:,d1)=sum(l,2)>0;
end
 l2{m}= l2{m}+1;
end

for m=1:size(Arcclass,2)
    
    for j=1:size(c,1)
        ven_data(j,m)=sum(l2{m}(:,1)==c(j,1) &...
                             l2{m}(:,2)==c(j,2) &...
                             l2{m}(:,3)==c(j,3))...
                             /size(l2{m},1);
    end
%     subplot(1,size(Arcclass,2)+1,m)
%     venn(ven_data(:,m));
%     axis equal off tight
%     hold all
end
    
    subplot(1,1,1)
    M=mean(ven_data,2);
    M
    [a Venn_data]=venn(M);
    axis equal off tight
    hold all
    for coor=1:7
        text(Venn_data.ZoneCentroid(coor,1),Venn_data.ZoneCentroid(coor,2)...
            ,num2str(round(M(coor,1)*100,1),3),'HorizontalAlignment', 'center','Fontsize',11,'FontName','Arial rounded mt bold')
    end
    
end
