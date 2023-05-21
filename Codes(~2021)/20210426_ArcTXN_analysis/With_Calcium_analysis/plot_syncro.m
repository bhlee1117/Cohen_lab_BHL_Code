function [d p_value]=plot_syncro(corr_mat,corr_mat_pooled,group_div,sprt,im_name,sw,rep_m,rep_d,xtick,g_show,t_mouse)
%% Show example adjacent matrix
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 300 300]);
if sw==1
    for t=1:size(corr_mat{rep_m,rep_d},3)
        ff=figure1;
        imagesc(corr_mat{rep_m,rep_d}(:,:,t),[0 0.3])
        axis equal tight
        colormap('jet')
        colorbar
        hold all
        for i=1:size(group_div{rep_m,rep_d},1)-1
            line([group_div{rep_m,rep_d}(i,1) group_div{rep_m,rep_d}(i,1)]+0.5,[0.5 size(corr_mat{rep_m,rep_d},1)+0.5],'color','w','linewidth',2)
            line([0.5 size(corr_mat{rep_m,rep_d},1)+0.5],[group_div{rep_m,rep_d}(i,1) group_div{rep_m,rep_d}(i,1)]+0.5,'color','w','linewidth',2)
        end
        F=getframe(ff);
        
        imwrite(F.cdata,im_name,'Writemode','Append')
    end
end
%%
if ~isempty(sprt)
    cond=[1 1;2 2]; cmap=distinguishable_colors(size(cond,1));
    figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
        'Color',[1 1 1],...
        'Renderer','painters','position',[100 100 400 250]);
    for cd=1:size(cond,1)
        lineProps.col{1}=cmap(cd,:);
        tmp=corr_mat_pooled{cond(cd,1),cond(cd,2)}{rep_m,rep_d};
        ave_corr=mean(reshape(tmp,size(tmp,1)*size(tmp,2),size(tmp,3)),1,'omitnan');
        std_corr=std(reshape(tmp,size(tmp,1)*size(tmp,2),size(tmp,3)),0,1,'omitnan')...
            ./sqrt(sum(~isnan(reshape(tmp,size(tmp,1)*size(tmp,2),size(tmp,3))),1));
        mseb([1:1:size(tmp,3)]*sprt,ave_corr,std_corr,lineProps,0.5)
        hold all
    end
    legend(xtick)
    xlabel('Time (sec)')
    ylabel('Average correlation')
else  % No segmentation
    
    %
    figure2 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
        'Color',[1 1 1],...
        'Renderer','painters','position',[100 100 300 300]);
    cmap=distinguishable_colors(size(corr_mat_pooled,1)*size(corr_mat_pooled,2));
    cmap2=[0.1 0.1 0.1;distinguishable_colors(3)];
    if size(g_show,1)==2
        G=[0.8 1.2];
    else if size(g_show,1)==3
            G=[0.7 1 1.3];
        else if size(g_show,1)==4
                G=[0.7 0.9 1.1 1.3];
            end
        end
    end
    
    for rep_d=1:3
        GG=G+rep_d-1;
    for i=1:size(g_show,1)
        if t_mouse
            for m=1:size(corr_mat_pooled{g_show(i,1),g_show(i,2)}(:,rep_d),1)
        d{rep_d,i}(m,1)=mean(corr_mat_pooled{g_show(i,1),g_show(i,2)}{m,rep_d},'omitnan');
            end
        else
        d{rep_d,i}=cell2mat(corr_mat_pooled{g_show(i,1),g_show(i,2)}(:,rep_d));
        end
        M(1,i)=mean(d{rep_d,i},'omitnan');
        S(1,i)=std(d{rep_d,i},0,1,'omitnan')/sqrt(sum(~isnan(d{rep_d,i})));
    end
%     boxplot(cell2mat(d(rep_d,:)),'positions', G+rep_d-1,'FullFactors','on')

%     violin(d(rep_d,:),'x',G+rep_d-1,'facecolor',cmap2,'edgecolor',[],'bw',0.1,'mc',[],'medc',[],'plotlegend',0);
    %

    %violin(d(rep_d,:),'facecolor',cmap,'edgecolor',[],'bw',0.1,'mc','[]','medc','[]','plotlegend',0);
     
%     b=bar(G+rep_d-1,M,'Barwidth',0.7,'LineWidth',2);
%     b.FaceColor='flat';
for j=1:size(g_show,1)
%b.CData(j,:)= cmap2(j,:);
plot(G(j)+rep_d-1,d{rep_d,j},'marker','o','color',cmap2(j,:),'linestyle','none')
    hold all
 end

    errorbar(G+rep_d-1,M,S,'LineWidth',2,'linestyle','none','color','k','Capsize',10,...
             'marker','+','markersize',10)
         hold all
         
    p_list=[1 2;2 3;1 3;3 4;2 4;1 4];
    if size(g_show,1)==2
        p_list=p_list(1,:); end
    if size(g_show,1)==3
        p_list=p_list(1:3,:); end
    if size(g_show,1)==4
        p_list=p_list(1:6,:); end
    for pp=1:size(p_list,1)
        if t_mouse
        [a p_value(pp,1)]=ttest(d{rep_d,p_list(pp,1)},d{rep_d,p_list(pp,2)});
         ref=1;
        else
        [a p_value(pp,1)]=kstest2(d{rep_d,p_list(pp,1)},d{rep_d,p_list(pp,2)});
        [a p_value(pp,2)]=ttest2(d{rep_d,p_list(pp,1)},d{rep_d,p_list(pp,2)});
        [p_value(pp,3) a]=ranksum(d{rep_d,p_list(pp,1)},d{rep_d,p_list(pp,2)});
         ref=3;
        end
    end
    p_value;
    axis tight
    set(gca,'XTick',[1:size(g_show,1)],'XTickLabel',xtick...
      ,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
    y_lim=0.05;
    for i=1:size(p_list,1)
        if p_value(i,ref)<0.05
            star='*';
            if p_value(i,ref)<0.01
                star='**';
                if p_value(i,ref)<0.001
                    star='***';
                end
            end
            line(GG(p_list(i,:)),[y_lim*0.72+y_lim*0.04*i y_lim*0.72+y_lim*0.04*i],'color','k','linewidth',2);
            line(GG([p_list(i,1) p_list(i,1)]),[y_lim*0.62+y_lim*0.04*i y_lim*0.72+y_lim*0.04*i],'color','k','linewidth',2);
            line(GG([p_list(i,2) p_list(i,2)]),[y_lim*0.62+y_lim*0.04*i y_lim*0.72+y_lim*0.04*i],'color','k','linewidth',2);
            text(mean(GG(p_list(i,:))),y_lim*0.74+y_lim*0.04*i,star,'FontSize',13,'FontName','arial rounded mt bold',...
                'HorizontalAlignment', 'center')
        end
    end
    xlim([0.3 rep_d+0.7])
    ylim([0 y_lim])
    ylabel('Correlation coefficient','FontName','arial rounded mt bold','FontSize',13)
    end
%         h = findobj(gca,'Tag','Box');
%         cmap2=repmat(cmap2,3,1);
% for j=1:length(h)
%     patch(get(h(j),'XData'),get(h(j),'YData'),cmap2(j,:),'FaceAlpha',.5);
% end

if size(g_show,1)==4
    set(gca,'Xtick',[repmat(G,1,3)+[0 0 0 0 1 1 1 1 2 2 2 2]],'Xticklabel',xtick)
else
    set(gca,'Xtick',[repmat(G,1,3)+[0 0 1 1 2 2]],'Xticklabel',xtick)
end
    xtickangle(45)
end

end