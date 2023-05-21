function P=plot_degrees(centrality,repd,t_mouse,sw,dot_sw,group_show,ticks)

% Draw some plots
% Degree centrality
if t_mouse
    if sw
        for d1=repd
            for g=group_show
                [ntmp NM]=cellfun(@size,centrality.dcdist{g}(:,d1),'UniformOutput',false);
                dN(:,g)=cell2mat(ntmp);
                dM{g}(:,d1)=cellfun(@mean,centrality.dcdist{g}(:,d1));         dM_n{g}(:,d1)=cellfun(@mean,centrality.dcnormdist{g}(:,d1));
                dS{g}(:,d1)=cellfun(@std,centrality.dcdist{g}(:,d1))./dN(:,g); dS_n{g}(:,d1)=cellfun(@std,centrality.dcnormdist{g}(:,d1))./dN(:,g);
            end
        end
          for g=1:size(centrality.dcdist,2)
            % plot_errorbar(dM{d1},'ttest',30,'Degrees',{'Arc^-','Arc^+'})
            plot_errorbar(dM_n{g},'ttest',1,'Normalized degrees',ticks,dot_sw)
        end
    else
        for g=group_show
            for d1=repd
                [ntmp NM]=cellfun(@size,centrality.dcdist{g}(:,d1),'UniformOutput',false);
                dN(:,g)=cell2mat(ntmp);
                dM{d1}(:,g)=cellfun(@mean,centrality.dcdist{g}(:,d1));         dM_n{d1}(:,g)=cellfun(@mean,centrality.dcnormdist{g}(:,d1));
                dS{d1}(:,g)=cellfun(@std,centrality.dcdist{g}(:,d1))./dN(:,g); dS_n{d1}(:,g)=cellfun(@std,centrality.dcnormdist{g}(:,d1))./dN(:,g);
            end
        end
        for d1=repd
            % plot_errorbar(dM{d1},'ttest',30,'Degrees',{'Arc^-','Arc^+'})
            plot_errorbar(dM_n{d1},'ttest',0.1,'Normalized degrees',ticks,dot_sw)
        end
    end
else
    if sw==1 % btw days
        cmap=distinguishable_colors(3);
        for g=group_show
            for d1=repd
                dc_pool{g,d1}=cell2mat(centrality.dcdist{g}(:,d1));
                dcnorm_pool{g,d1}=cell2mat(centrality.dcnormdist{g}(:,d1));
            end
            P=plot_errorbar(dcnorm_pool(g,:),'ranksum',0.2,'Normalized degrees',ticks,dot_sw)
            
        end
    else     % btw the groups
        figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],'Renderer','painters','position',[100 100 500 300]);
        cmap=[0 0 0;1 0 0];
        if length(group_show)==2
            grd=[0.8 1.2];
        else if length(group_show)==3
                grd=[0.7 1 1.3];
            else if length(group_show)==4
                    grd=[0.7 0.9 1.1 1.3];
                end
            end
        end
        GG=[];
        for d1=repd
            for g=group_show
                dc_pool{g,d1}=cell2mat(centrality.dcdist{g}(:,d1));
                dcnorm_pool{g,d1}=cell2mat(centrality.dcnormdist{g}(:,d1));
            end
            plot_errorbar2(grd+d1-1,dcnorm_pool(:,d1)','ttest2',0.3,'Normalized degrees',ticks,dot_sw)
            xtickangle(45)
            GG=[GG grd+d1-1];
        end
        
        set(gca,'XTick',GG,'XTickLabel',ticks...
    ,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
        xlim([0.3 d1+0.7])
    end
end
% %%
% % % Degree histograms
% bin=[0:0.02:0.2];
% cmap=distinguishable_colors(3);
% clear dc_histo
% if t_mouse
%     if sw==1 % btw days
%         for g=1:size(centrality.dcdist,2)
%             figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
%                 'Color',[1 1 1],'Renderer','painters','position',[100 100 300 300]);
%             for d1=repd
%                 for m=1:size(centrality.dcnormdist{g},1)
%                     figure(1)
%                     h=histogram(centrality.dcnormdist{g}{m,d1},bin,'Normalization','probability');
%                     dc_histo{g,d1}(m,:)=h.Values;
%                     close(figure(1))
%                 end
%                 lineProps.col{1}=cmap(d1,:);
%                 mseb(bin(1:end-1),mean(dc_histo{g,d1},1,'omitnan'),std(dc_histo{g,d1},0,1,'omitnan')./sqrt(sum(~isnan(dc_histo{g,d1}),1)),lineProps,1)
%                 hold all
%             end
%             ylim([0 0.7])
%             ylabel('Probability','FontName','arial rounded mt bold','FontSize',13)
%             xlabel('Normalized degree','FontName','arial rounded mt bold','FontSize',13)
%             set(gca,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
%             legend({'Day 1','Day 2','Day 3'});
%         end
%         
%     else  % btw the groups
%         for d1=repd
%             figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
%                 'Color',[1 1 1],'Renderer','painters','position',[100 100 300 300]);
%             for g=1:size(centrality.dcdist,2)
%                 for m=1:size(centrality.dcnormdist{g},1)
%                     h=histogram(centrality.dcnormdist{g}{m,d1},bin,'Normalization','probability');
%                     dc_histo{g,d1}(m,:)=h.Values;
%                 end
%                 lineProps.col{1}=cmap(d1,:);
%                 mseb(bin,mean(dc_histo{g,d1},1,'omitnan'),std(dc_histo{g,d1},0,1,'omitnan')./sqrt(sum(isnan(dc_histo{g,d1}),1)),lineProps,1)
%                 hold all
%             end
%         end
%     end
% else
%     if sw
%     else
%     end
% end



end
