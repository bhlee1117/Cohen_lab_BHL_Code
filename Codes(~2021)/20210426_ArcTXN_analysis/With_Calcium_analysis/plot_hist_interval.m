function plot_hist_interval(I,day,bins,ticks,y_lim,mthd)
for d1=day
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 450 350]);
M=[]; S=[];
for g=1:size(I,2)

histogram(cell2mat(I{g}(:,d1)),bins,'Normalization',mthd,'linewidth',0.5)
hold all
% if t_mouse
% M=[M mean(cell2mat(I{g}(:,d1)),1,'omitnan')];
% S=[S  std(cell2mat(I{g}(:,d1)),0,1,'omitnan')/sqrt(size(cell2mat(I{g}(:,d1)),1))];    
% else
M=[M mean(cell2mat(I{g}(:,d1)),1,'omitnan')];
S=[S  std(cell2mat(I{g}(:,d1)),0,1,'omitnan')/sqrt(size(cell2mat(I{g}(:,d1)),1))];
% end
end
p_list=[1 2;1 3;2 3;1 4;2 4;3 4];
if size(I,2)==2
    p_list=p_list(1,:); end
    if size(I,2)==3
    p_list=p_list(1:3,:); end
    if size(I,2)==4
    p_list=p_list(1:6,:); end
for i=1:size(p_list,1) 
[p_value(i,1) a]=ranksum(cell2mat(I{p_list(i,1)}(:,d1)),cell2mat(I{p_list(i,2)}(:,d1)));
[a p_value(i,2)]=kstest2(cell2mat(I{p_list(i,1)}(:,d1)),cell2mat(I{p_list(i,2)}(:,d1)));
[a p_value(i,3)]=ttest2(cell2mat(I{p_list(i,1)}(:,d1)),cell2mat(I{p_list(i,2)}(:,d1)));
end
p_value
set(gca,'FontSize',8,'LineWidth',1,...
   'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
ylabel('Probability','LineWidth',2,'FontSize',13,...
    'FontName','arial rounded mt bold')
xlabel('Ca^2^+ interval (s)','LineWidth',2,'FontSize',13,...
    'FontName','arial rounded mt bold')
legend(ticks)

figure2 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 300 300]);

errorbar([1:size(I,2)],M,S,...
            'LineWidth',2,'linestyle','-','color','k','Capsize',10)
        hold all
        set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:size(I,2)],'XTickLabel',...
            ticks,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
        ylabel('Ca^2^+ interval (s)','LineWidth',2,'FontSize',13,...
            'FontName','arial rounded mt bold')
        xlim([0.3 size(I,2)+.7])
        ylim([0 y_lim])
        ref=1;
                for i=1:size(p_list,1)
            if p_value(i,ref)<0.05
                star='*';
                if p_value(i,ref)<0.01
                    star='**';
                    if p_value(i,ref)<0.001
                        star='***';
                    end
                end
                
                line(p_list(i,:),[y_lim*0.72+y_lim*0.04*i y_lim*0.72+y_lim*0.04*i],'color','k','linewidth',2);
                line([p_list(i,1) p_list(i,1)],[y_lim*0.62+y_lim*0.04*i y_lim*0.72+y_lim*0.04*i],'color','k','linewidth',2);
                line([p_list(i,2) p_list(i,2)],[y_lim*0.62+y_lim*0.04*i y_lim*0.72+y_lim*0.04*i],'color','k','linewidth',2);
                text(mean(p_list(i,:)),y_lim*0.74+y_lim*0.04*i,star,'FontSize',13,'FontName','arial rounded mt bold',...
                    'HorizontalAlignment', 'center')
            end
        end
        xtickangle(45)
end
end