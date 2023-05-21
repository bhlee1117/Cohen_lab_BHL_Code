%function plot_estrada(estrada_index,repd,t_mouse,sw,dot_sw)
%%
% Draw some plots
% Degree centrality
if t_mouse
    for g=1:size(estrada_index,2)
        for d1=repd
            for m=1:size(estrada_index{g},1)
            dM{d1}(m,g)=mean(estrada_index{g}{m,d1},'omitnan');         
            dS{d1}(m,g)=std(estrada_index{g}{m,d1},0,1,'omitnan')./sum(~isnan(estrada_index{g}{m,d1}));
            end
        end
    end
    for d1=repd
        plot_errorbar(dM{d1},'ttest',0.07,'Normalized SC',{'Arc^-','Arc^+'},dot_sw)
    end
else
    if sw==1 % in each day
        for g=1:size(estrada_index,2)
            for d1=repd
                ee_pool{g,d1}=cell2mat(estrada_index{g}(:,d1));
            end
            plot_errorbar(ee_pool(g,:),'ttest2',0.07,'Normalized SC',{'Day 1','Day 2','Day 3'},dot_sw)
        end
    else     % btw the groups
        for d1=repd
            for g=1:size(estrada_index,2)
                ee_pool{g,d1}=cell2mat(estrada_index{g}(:,d1));
            end
            plot_errorbar(ee_pool(:,d1)','ttest2',0.07,'Normalized SC',{'Arc^-','Arc^+'},dot_sw)
        end
    end
end
%
% % Degree histograms
% for
% end



% end
