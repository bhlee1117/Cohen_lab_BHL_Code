function plot_dff_astro(C_df,timescale,stim_time,ticks,color)
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 400 300]);
ylim=5;
p=[];
for i=1:size(C_df,1) %group
    for j=1:size(C_df,2) % mouse
        lineProps.col{1}=color(size(C_df,2)*(i-1)+j,:);
        pp=mseb([0:timescale:timescale*(size(C_df{i,j},2)-1)],mean(C_df{i,j},1),std(C_df{i,j},0,1)/sqrt(size(C_df{i,j},1)),lineProps,0.5);
        p=[p pp.mainLine];
        hold all
    end
end

fill([stim_time  timescale*(size(C_df{i,j},2)-1)  timescale*(size(C_df{i,j},2)-1) stim_time], [0 0 ylim ylim], 'r','Linestyle','none');
alpha(0.2)
legend(p, ticks,'Location','northwest','FontName','arial rounded mt bold','FontSize',10,'LineWidth',1)
ylabel('\Delta F/F','FontSize',13,'FontName','arial rounded mt bold')
xlabel('Time (sec)','FontSize',13,'FontName','arial rounded mt bold')

end