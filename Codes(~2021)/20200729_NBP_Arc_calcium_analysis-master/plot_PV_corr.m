function p=plot_PV_corr(PV_corr,xtick,sw,t_mouse)
cmap=[0.2 0.2 0.2;1 0 0];
if sw==0
    figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
        'Color',[1 1 1],...
        'Renderer','painters','position',[100 100 300 300]);
    M=mean(PV_corr,'omitnan');
    S=std(PV_corr,0,1,'omitnan')./sqrt(sum(~isnan(PV_corr),1));
    pp=[1 2;1 3;2 3];
    if size(PV_corr,2)<3
        pp=pp(1,:);
    end
    for i=1:size(pp,1)
        [a p(i)]=ttest(PV_corr(:,pp(i,1)),PV_corr(:,pp(i,2)));
    end
    errorbar([1:size(PV_corr,2)],M,S, 'LineWidth',2,'linestyle','none','color','k','Capsize',10,'marker','+','markersize',10)
    hold all
    % b=bar([1 2 3],M,'Barwidth',0.7,'LineWidth',2);
    % b.FaceColor='flat';
    % for i=1:3
    % b.CData(i,:)= cmap(i,:);
    % hold all
    % %plot(X,dat{i},'marker','.','markersize',10,'linestyle','none','color',cmap(i,:)/2)
    % end
    set(gca,'FontSize',8,'LineWidth',1,'XTick',[1 2 3],'XTickLabel',...
        xtick,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
    ylabel('Spatial correlation','LineWidth',2,'FontSize',13,...
        'FontName','arial rounded mt bold')
    xlim([0.5 size(PV_corr,2)+.5])
    ylim([-0.1 0.55])
    for i=1:size(pp,1)
        if p(i)<0.05
            star='*';
            if p(i)<0.01
                star='**';
                if p(i)<0.001
                    star='***';
                end
            end
            
            line(pp(i,:),[0.4+0.04*i 0.4+0.04*i],'color','k','linewidth',2);
            line([pp(i,1) pp(i,1)],[0.34+0.04*i 0.4+0.04*i],'color','k','linewidth',2);
            line([pp(i,2) pp(i,2)],[0.34+0.04*i 0.4+0.04*i],'color','k','linewidth',2);
            text(mean(pp(i,:)),0.41+0.04*i,star,'FontSize',13,'FontName','arial rounded mt bold',...
                'HorizontalAlignment', 'center')
        end
    end
else
    if t_mouse
        d=cellfun(@mean,PV_corr);
        M=mean(d,1); S=std(d,0,1)/sqrt(size(d,1));
        [a p]=ttest(d(:,1),d(:,2))
        plot_errorbar(d,'ttest',0.7,'Spatial correlation',{'Arc^-','Arc^+'},0)
    hold all
    b=bar([1 2],M,'Barwidth',0.7,'LineWidth',2);
    b.FaceColor='flat';
    else
        M=cellfun(@mean,PV_corr);
        S=[std(PV_corr{1})/sqrt(size(PV_corr{1},1)) std(PV_corr{2})/sqrt(size(PV_corr{2},1))];
        [a p(1)]=kstest2(PV_corr{1},PV_corr{2});
        [a p(2)]=ttest2(PV_corr{1},PV_corr{2});
        [p(3) a]=ranksum(PV_corr{1},PV_corr{2});
        plot_errorbar(PV_corr,'ranksum',0.7,'Spatial correlation',{'Arc^-','Arc^+'},0)
    hold all
    b=bar([1 2],M,'Barwidth',0.7,'LineWidth',2);
    b.FaceColor='flat';
    end
    
    for i=1:2
        b.CData(i,:)= cmap(i,:);
    end
end
end
