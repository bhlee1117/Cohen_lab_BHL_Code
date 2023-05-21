function p_value=plot_errorbar2(GG,dddd,method,y_lim,ylab,xtick,dot_sw)
% figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
%     'Color',[1 1 1],'Renderer','painters','position',[100 100 300 300]);
cmap=[0.1 0.1 0.1; 1 0 0; distinguishable_colors(5)];
p_list=[1 2;2 3;1 3;3 4;2 4;1 4];
if size(dddd,2)==2
    p_list=p_list(1,:); end
if size(dddd,2)==3
    p_list=p_list(1:3,:); end
if size(dddd,2)==4
    p_list=p_list(1:6,:); end

if iscell(dddd)
    M=[]; S=[];
    for g=1:size(dddd,2)
        M=[M mean(dddd{g},1,'omitnan')];
        S=[S  std(dddd{g},'omitnan')/sqrt(size(dddd{g},1))];
    end
%     b=bar(GG,M,'Barwidth',0.7,'LineWidth',2);
     hold all
%     b.FaceColor='flat';
%     for j=1:size(dddd,2)
%         b.CData(j,:)= cmap(j,:);
%     end
    for g=1:size(dddd,2)
        G=GG(g)+randn(size(dddd{g},1),1)*0.02;
        if dot_sw
            for u=1:size(dddd{g},1)
                plot(G,dddd{g},'marker','.','color',[0.7 0.7 0.7],'markersize',5,'linestyle','none')
                hold all
            end
        end
    end
    
    errorbar(GG,M,S,'LineWidth',2,'linestyle','none',...
        'color','k','Capsize',10,'marker','+','markersize',10)
    xlim([GG(1,1)-0.2 GG(end)+0.2])
else
    if dot_sw
        for i=1:size(dddd,1)
            plot(GG,dddd(i,:),'color',[0.7 0.7 0.7],'LineWidth',1.5,'marker','o')
            hold all
        end
    end
    errorbar(GG,mean(dddd,1,'omitnan'),std(dddd,0,1,'omitnan')./sqrt(sum(~isnan(dddd),1)),'LineWidth',2,'linestyle','none',...
        'color','k','Capsize',10,'marker','+','markersize',10)
    hold all
    
%     b=bar(GG,mean(dddd,1,'omitnan'),'Barwidth',0.7,'LineWidth',2);
%     hold all
%     b.FaceColor='flat';
%     for j=1:size(dddd,2)
%         b.CData(j,:)= cmap(j,:);
%     end
    xlim([GG(1,1)-0.2 GG(end)+0.2])
end
for pp=1:size(p_list,1)
    switch method
        case 'ttest'
            [a p_value(pp,1)]=ttest(dddd(:,p_list(pp,1)),dddd(:,p_list(pp,2)));
        case 'ttest2'
            [a p_value(pp,1)]=ttest2(dddd{p_list(pp,1)},dddd{p_list(pp,2)});
        case 'kstest2'
            [a p_value(pp,1)]=kstest2(dddd{p_list(pp,1)},dddd{p_list(pp,2)});
        case 'ranksum'
            [p_value(pp,1) a]=ranksum(dddd{p_list(pp,1)},dddd{p_list(pp,2)});
    end
end
ref=1;
i=1;
for i=1:size(p_list,1)
    if p_value(i,ref)<0.05
        star='*';
        if p_value(i,ref)<0.01
            star='**';
            if p_value(i,ref)<0.001
                star='***';
            end
        end
        line(GG([p_list(i,:)]),[y_lim*0.72+y_lim*0.04*i y_lim*0.72+y_lim*0.04*i],'color','k','linewidth',2);
        line(GG([p_list(i,1) p_list(i,1)]),[y_lim*0.62+y_lim*0.04*i y_lim*0.72+y_lim*0.04*i],'color','k','linewidth',2);
        line(GG([p_list(i,2) p_list(i,2)]),[y_lim*0.62+y_lim*0.04*i y_lim*0.72+y_lim*0.04*i],'color','k','linewidth',2);
        text(mean(GG(p_list(i,:))),y_lim*0.74+y_lim*0.04*i,star,'FontSize',13,'FontName','arial rounded mt bold',...
            'HorizontalAlignment', 'center')
    end
end
ylim([0 y_lim])
ylabel(ylab,'FontName','arial rounded mt bold','FontSize',13)
set(gca,'XTick',GG,'XTickLabel',xtick...
    ,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
end