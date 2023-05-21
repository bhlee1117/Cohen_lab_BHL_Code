function fourier=Plot_fft_cond(dat,groups,xtick,y_lim,x_lim,sw,post_ref,t_mouse)


for m=1:size(dat,2)
    
    for d1=1:3
        for g=1:size(groups,2)
            fourier{g}{m,d1}=[];
        end
    end
    
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
            fourier{g}{m,d1}=[];
            for l=find(list{m,g})'
                if size(dat{m}.Cal{l,d1},2)>1 
                  fourier{g}{m,d1}=[fourier{g}{m,d1}; dat{m}.FT{l,3+d1}];
                else
                  fourier{g}{m,d1}=[fourier{g}{m,d1}; NaN(1,6001)]; 
                end
            end
        end
    end
end

%% plot

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
                    MF(m,:)=mean(fourier{g}(m,d1),1,'omitnan'); SF(m,:)=mean(fourier{g}(m,d1),0,1,'omitnan'); end
                mseb([0:15/6000:15],mean(MF,1,'omitnan'),std(SF,0,1,'omitnan')/sqrt(size(SF,1)),lineProps,0.5);
            else
                datum=cell2mat(fourier{g}(:,d1));
                mseb([0:15/6000:15],mean(datum,1,'omitnan'),std(datum,0,1,'omitnan')/sqrt(size(datum,1)),lineProps,0.5);
            end
        end
        hold all
        
    end
    ylabel('Auto-correlation')
    xlabel('Time (sec)')
    xlim([0 x_lim])
    ylim([0 y_lim])
    legend(xtick)
    
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
                    MF(m,:)=mean(fourier{g}(m,d1),1,'omitnan'); SF(m,:)=mean(fourier{g}(m,d1),0,1,'omitnan'); end
                mseb([0:15/6000:15],mean(MF,1,'omitnan'),std(MF,0,1,'omitnan')/sqrt(size(MF,1)),lineProps,0.5);
            else
                datum=cell2mat(fourier{g}(:,d1));
                mseb([0:15/6000:15],mean(datum,1,'omitnan'),std(datum,0,1,'omitnan')/sqrt(size(datum,1)),lineProps,0.5);
            end
            hold all
        end
        ylabel('Power spectrum')
        xlabel('Frequency (sec)')
        xlim([0 x_lim])
        ylim([0.001 y_lim])
        legend(xtick)
        set(gca,'yscale','log','FontName','arial rounded mt bold','FontSize',13,'LineWidth',2)
    end
end
end