 function plot_graphs(corr_mat,group_div,corr_th,repd)
%
cmap_group=[0 0 0;1 0 0];
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 400 400]);

        A=round(corr_mat{repd(1,1),repd(1,2)},3);
        A(find(diag(zeros(size(A,1),1)+1)==1))=0;
        A=double(A>corr_th);
        G=graph(double(A));
        dc_norm=sum(A,2)/sum(A(:))*2;
            gplot=plot(G,'markersize',5,'EdgeColor',[0 0 0],'NodeLabel',{});
            hold all
            cmap_g=jet(round((max(dc_norm)-min(dc_norm))*1000)+1);
            gplot.NodeColor=cmap_g(round((dc_norm-min(dc_norm))*1000)+1,:);
            divide=[[1;group_div{repd(1,1),repd(1,2)}(1:end-1,1)+1] [group_div{repd(1,1),repd(1,2)}]];
            for g=1:size(group_div{repd(1,1),repd(1,2)},1)
                plot(gplot.XData(divide(g,1):divide(g,2)),gplot.YData(divide(g,1):divide(g,2)),'o','markersize',6,'linewidth',1.5,'color',cmap_group(g,:))
                axis equal
            end
     
 end