function plot_cluster_coeff(cc_dist,repd,t_mouse,sw,dot_sw,group_show,ticks)
%%
% Draw some plots
% Degree centrality
if t_mouse
    for g=group_show
        for d1=repd
            for m=1:size(cc_dist{g},1)
            dM{d1}(m,g)=mean(cc_dist{g}{m,d1},'omitnan');         
            dS{d1}(m,g)=std(cc_dist{g}{m,d1},0,1,'omitnan')./sum(~isnan(cc_dist{g}{m,d1}));
            end
        end
    end
    for d1=repd
        plot_errorbar(dM{d1},'ttest',1,'Cluster coefficient',ticks,dot_sw)
    end
else
    if sw==1 % in each day
        for g=group_show
            for d1=repd
                cc_pool{g,d1}=cell2mat(cc_dist{g}(:,d1));
            end
            plot_errorbar(cc_pool(g,:),'ranksum',0.6,'Cluster coefficient',ticks,dot_sw)
        end
    else     % btw the groups
        figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],'Renderer','painters','position',[100 100 500 300]);
        if length(group_show)==2
            grd=[0.8 1.2];
        else if length(group_show)==3
                grd=[0.7 1 1.3];
            else if length(group_show)==4
                    grd=[0.7 0.9 1.1 1.3];
                end
            end
        end
        cmap=[0 0 0;1 0 0];
        GG=[];
        for d1=repd
            for g=group_show
                cc_pool{g,d1}=cell2mat(cc_dist{g}(:,d1));
            end
            plot_errorbar2(grd+d1-1,cc_pool(:,d1)','kstest2',1,'Cluster coefficient',ticks,dot_sw)
            xtickangle(45)
             GG=[GG grd+d1-1];
        end
        set(gca,'XTick',GG,'XTickLabel',ticks...
    ,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
        xlim([0.3 d1+0.7])
    end
end

%
% % Degree histograms
% for
% end



 end
