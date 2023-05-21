function [IBI_cond datum p_value]=burst_detection_cond(IBI,dat,groups,sw,x_lim,y_lim,t_mouse,xtick,post_ref)
void_th=0;
D_th=2.5;
method='ttest2';
clear list IBI_cond
for m=1:size(dat,2)
        for d1=1:3
            IBI.Brst{m,d1}=[IBI.Brst{m,d1} IBI.theta_Brst{m,d1}];
        end
end
for m=1:size(dat,2)
    for g=1:size(groups,2)
        for d1=1:3
            IBI_cond{g,d1}{m,1}=[];
            Pow_cond{g,d1}{m,1}=[];
            Brst_cond{g,d1}{m,1}=[];
            theta_Brst_cond{g,d1}{m,1}=[]; end
    end
end

for m=1:size(dat,2)
    if post_ref
        Arc_class{m}=[dat{m}.Arc_post_ref(:,2:end)];
    else
        Arc_class{m}=[dat{m}.Arc(:,2:end)];
    end
    
    for g=1:size(groups,2)
        for cd=1:size(groups{g},2)/2
            list{m,g}(:,cd)=Arc_class{m}(:,groups{g}(1,2*cd-1))==groups{g}(1,2*cd);
        end
        list{m,g}=sum(list{m,g},2)==size(groups{g},2)/2; % cells that match the conditions
        
        for d1=1:3
            IBI_cond{g,d1}{m,1}=[IBI_cond{g,d1}{m,1}; IBI.IBI{m,d1}(find(list{m,g}),:)];
            Pow_cond{g,d1}{m,1}=[Pow_cond{g,d1}{m,1}; IBI.Power_spec{m,d1}(find(list{m,g}),:)];
            Brst_cond{g,d1}{m,1}=[Brst_cond{g,d1}{m,1}; IBI.Brst{m,d1}(find(list{m,g}),:)];
            theta_Brst_cond{g,d1}{m,1}=[theta_Brst_cond{g,d1}{m,1}; IBI.theta_Brst{m,d1}(find(list{m,g}),:)];
        end
    end
end

%% Now plot
% Inter-burst interval
if sw
    cmap=distinguishable_colors(3);
    for g=1:size(groups,2)
        figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
            'Color',[1 1 1],...
            'Renderer','painters','position',[100 100 400 250]);
        
        for d1=1:3
            lineProps.col{1}=cmap(d1,:);
            if t_mouse
                for m=1:size(dat,2)
                MF(m,:)=mean(IBI_cond{g,d1}{m},1,'omitnan'); end
                mseb([1:size(MF,2)]/30,mean(MF,1,'omitnan'),std(MF,0,1,'omitnan')/sqrt(size(MF,1)),lineProps,0.5);
            else
                datum=cell2mat(IBI_cond{g,d1});
                mseb([1:size(datum,2)]/30,mean(datum,1,'omitnan'),std(datum,0,1,'omitnan')./sqrt(sum(~isnan(datum),1)),lineProps,1);
            end
        end
        hold all
        ylabel('Burst count','FontName','arial rounded mt bold','FontSize',13)
        xlabel('Inter-burst interval (s)','FontName','arial rounded mt bold','FontSize',13)
        xlim([0 x_lim])
        ylim([0 y_lim])
        legend(xtick)
    end
    
    
else
    cmap=distinguishable_colors(size(groups,2));
    for d1=1:3
        figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
            'Color',[1 1 1],...
            'Renderer','painters','position',[100 100 400 250]);
        
        for g=1:size(groups,2)
            lineProps.col{1}=cmap(g,:);
            if t_mouse
                for m=1:size(dat,2)
                    MF(m,:)=mean(IBI_cond{g,d1}{m},1,'omitnan');  end
                mseb([1:size(MF,2)]/30,mean(MF,1,'omitnan'),std(MF,0,1,'omitnan')/sqrt(size(MF,1)),lineProps,0.5);
            else
                datum=cell2mat(IBI_cond{g,d1});
                mseb([1:size(datum,2)]/30,mean(datum,1,'omitnan'),std(datum,0,1,'omitnan')./sqrt(sum(~isnan(datum),1)),lineProps,1);
            end
            hold all
        end
        ylabel('Burst probability','FontName','arial rounded mt bold','FontSize',13)
        xlabel('Inter-burst interval (s)','FontName','arial rounded mt bold','FontSize',13)
        xlim([0 x_lim])
        ylim([0 y_lim])
        legend(xtick)
        legend(xtick)
        
    end
end

% Power spectrum
clear MF
if sw
    cmap=distinguishable_colors(3);
    for g=1:size(groups,2)
        figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
            'Color',[1 1 1],...
            'Renderer','painters','position',[100 100 400 250]);
        
        for d1=1:3
            lineProps.col{1}=cmap(d1,:);
            if t_mouse
                for m=1:size(dat,2)
                MF(m,:)=mean(Pow_cond{g,d1}{m},1,'omitnan');  end
                mseb([0:IBI.f(2):IBI.f(2)*(size(MF,2)-1)],mean(MF,1,'omitnan'),std(MF,0,1,'omitnan')/sqrt(size(MF,1)),lineProps,0.5);
            else
                datum=cell2mat(Pow_cond{g,d1});
                mseb([0:IBI.f(2):IBI.f(2)*(size(datum,2)-1)],mean(datum,1,'omitnan'),std(datum,0,1,'omitnan')./sqrt(sum(~isnan(datum),1)),lineProps,1);
            end
        end
        hold all
        ylabel('Burst probability','FontName','arial rounded mt bold','FontSize',13)
        xlabel('Inter-burst interval (s)','FontName','arial rounded mt bold','FontSize',13)
        %     xlim([0 x_lim])
        %     ylim([0 y_lim])
        legend(xtick)
        set(gca,'yscale','log')
    end
    
    
else
    cmap=distinguishable_colors(size(groups,2));
    for d1=1:3
        figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
            'Color',[1 1 1],...
            'Renderer','painters','position',[100 100 400 250]);
        
        for g=1:size(groups,2)
            lineProps.col{1}=cmap(g,:);
            if t_mouse
                clear MF
                for m=1:size(dat,2)
                    MF(m,:)=mean(Pow_cond{g,d1}{m},1,'omitnan');  end
                mseb([0:IBI.f(2):IBI.f(2)*(size(MF,2)-1)],mean(MF,1,'omitnan'),std(MF,0,1,'omitnan')/sqrt(size(MF,1)),lineProps,0.5);
            else
                datum=cell2mat(Pow_cond{g,d1});
                mseb([0:IBI.f(2):IBI.f(2)*(size(datum,2)-1)],mean(datum,1,'omitnan'),std(datum,0,1,'omitnan')./sqrt(sum(~isnan(datum),1)),lineProps,1);
            end
            hold all
        end
        ylabel('Power (A.U.)','FontName','arial rounded mt bold','FontSize',13)
        xlabel('Frequency (Hz)','FontName','arial rounded mt bold','FontSize',13)
        %         xlim([0 x_lim])
        %         ylim([0 y_lim])
        legend(xtick)
        set(gca,'yscale','log')
    end
end


% Burst properties
ylabs={'Burst rate (Hz)','Mean spike number (/burst)','Max spike number (/burst)','Theta burst rate (s^{-1})'};
clear MF
if sw
    cmap=distinguishable_colors(3);
    G=[0.7 1 1.3];
    for cat=4
          GG=[];  yl=0;
        figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
            'Color',[1 1 1],...
            'Renderer','painters','position',[100 100 400 250]);
        if t_mouse
            for g=1:size(groups,2)
                for d1=1:3
                    for m=1:size(dat,2)
                   MF{cat,1}(m,d1)=mean(Brst_cond{g,d1}{m}(:,cat),1,'omitnan'); end
                end
                GG=[GG G+g-1];
                 if max(cell2mat(cellfun(@mean_omitnan,MF(cat,:),'UniformOutput',0)))*1.5>yl
                    yl=max(cell2mat(cellfun(@mean_omitnan,MF(cat,:),'UniformOutput',0)))*1.5;
                end
                plot_errorbar2(G+g-1,MF{cat},'ttest',yl,ylabs{cat},xtick,1)
            end
            set(gca,'Xtick',GG,'XtickLabel',xtick)
            
        else
            clear datum
            for g=1:size(groups,2)
                for d1=1:3
                    tmp=cell2mat(Brst_cond{g,d1});
                    datum{cat}{d1}=tmp(:,cat);
                end
                GG=[GG G+1-1];
                 if max(cellfun(@mean_omitnan,datum{cat}))*1.5>yl
                    yl=max(cellfun(@mean_omitnan,datum{cat}))*1.5;
                end
                p_value=plot_errorbar2(G+1-1,datum{cat},method,yl,ylabs{cat},xtick,0)
            end
            set(gca,'Xtick',G,'XtickLabel',xtick)
            xlim([GG(1,1)-0.2 GG(end)+0.2])
            xtickangle(45)
        end
    end
else
    if size(groups,2)==2
        G=[0.8 1.2];
    else if size(groups,2)==3
            G=[0.7 1 1.3];
        else
            G=[0.7 0.9 1.1 1.3];
        end
    end
    cmap=distinguishable_colors(size(groups,2)); clear MF
    for cat=4
        GG=[]; yl=0;
        figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
            'Color',[1 1 1],...
            'Renderer','painters','position',[100 100 400 250]);
        for d1=1:3
            if t_mouse
                for g=1:size(groups,2)
                    for m=1:size(dat,2)
                        MF{cat}(m,g)=mean(Brst_cond{g,d1}{m}(:,cat),1,'omitnan'); end
                end
                GG=[GG G+d1-1];
                if max(mean(MF{cat}))*2>yl
                    yl=max(mean(MF{cat}))*2;
                end
                plot_errorbar2(G+d1-1,MF{cat},'ttest',yl,ylabs{cat},xtick,1)
                hold all
            else
                clear datum
                for g=1:size(groups,2)
                    tmp{g}=cell2mat(Brst_cond{g,d1});
                    datum{g}=tmp{g}(:,cat);
                end
                  if max(cellfun(@mean_omitnan,datum))*1.5>yl
                    yl=max(cellfun(@mean_omitnan,datum))*1.5;
                end
                GG=[GG G+d1-1];
                p_value=plot_errorbar2(G+d1-1,datum,method,yl,ylabs{cat},xtick,0)
                
            end
        end
        set(gca,'Xtick',GG,'XtickLabel',xtick)
        xlim([GG(1,1)-0.2 GG(end)+0.2])
        xtickangle(45)
    end
end


end
function M=mean_omitnan(d)
M=mean(d,1,'omitnan');
end
