function cal_im_show(DFF,group,mouse,timescale)
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 500 300]);
imagesc(DFF{group,mouse},[0 3])
colormap('jet')
axis tight equal
set(gca,'XTick',[20:20:size(DFF{group,mouse},2)],'XTickLabel',round([timescale*20:timescale*20:timescale*size(DFF{group,mouse},2)]),'FontName','arial rounded mt bold','FontSize',10,'LineWidth',2,...
              'YTick',[8:8:size(DFF{group,mouse},1)])
xlabel('Time (sec)')
ylabel('Cell #')
colorbar
end