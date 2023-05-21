function [HR]=plot_histo_bin_trace_nogroup(dat,tick,bin,t_mouse)

for m=1:size(dat,2) %mouse
    for cd=1:4
        HR{m,cd}=[];
    end
    
    ylab='Ca^2^+ event rate (Hz)';
    for d1=1:3
        [r c]=cellfun(@size,dat{m}.Cal(:,3+d1));  max_bin=ceil(max(c)*1/30/bin);
        for c=1:size(dat{m}.Peak,1)
            if size(dat{m}.Peak{c,d1},2)==1
                binned_tr{m}{c,d1}=NaN;
            else
                %binned_tr{m}{c,d1}=peak_to_trace(dat{m}.Peak{c,d1},dat{m}.Cal{c,3+d1},bin,max_bin);
                binned_tr{m}{c,d1}=movmean(dat{m}.S{c,d1},30,'omitnan')*30;
                %binned_tr{m}{c,d1}=dat{m}.S{c,3+d1};
                binned_tr{m}{c,d1}(binned_tr{m}{c,d1}==0)=NaN;
            end
        end
    end
    datum=binned_tr{m}(:);
    %     if post_ref
    %         Arc_class{m}=[dat{m}.Arc_post_ref(:,2:end)];
    %     else
    Arc_class{m}=[dat{m}.Arc(:,2:end)];
    Arc_class{m}=Arc_class{m}(:);
    %     end
    for cd=1:4
        list{m,cd}=Arc_class{m}==cd;
        for c1=find(list{m,cd})'
            if size(datum{c1},2)==1
                h=NaN(1,30);
            else
                h=histcounts(datum{c1},[0:1:30],'Normalization','Probability');
            end
            HR{m,cd}=[HR{m,cd};h];
        end
    end
end
%% Now plot
if t_mouse
    cmap=distinguishable_colors(4);
    figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
        'Color',[1 1 1],'Renderer','painters','position',[100 100 300 200]);
    
    for cd=1:4
        for m=1:size(dat,2)
            HM(m,:)=mean(HR{m,cd},1,'omitnan');
        end
        lineProps.col{1}=cmap(cd,:);
        mseb([0:size(HM,2)-1]*bin,mean(HM,1,'omitnan')...
            ,std(HM,0,1,'omitnan')/sqrt(size(HM,1)),lineProps,1);
        hold all
    end
    legend(tick)
    xlabel(ylab,'LineWidth',2,'FontSize',13,...
        'FontName','arial rounded mt bold')
    ylabel('Probability','LineWidth',2,'FontSize',13,...
        'FontName','arial rounded mt bold')
    set(gca,'xscale','linear','yscale','log','FontSize',8,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
else
    cmap=distinguishable_colors(4);
    
    
    figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
        'Color',[1 1 1],'Renderer','painters','position',[100 100 300 200]);
    
    for cd=1:4
        Htmp=cell2mat(HR(:,cd));
        lineProps.col{1}=cmap(cd,:);
        mseb([0:size(Htmp,2)-1]*bin,mean(Htmp,1,'omitnan')...
            ,std(Htmp,0,1,'omitnan')./sqrt(sum(~isnan(Htmp),1)),lineProps,1);
        hold all
    end
    legend(tick)
    xlabel(ylab,'LineWidth',2,'FontSize',13,...
        'FontName','arial rounded mt bold')
    ylabel('Probability','LineWidth',2,'FontSize',13,...
        'FontName','arial rounded mt bold')
    set(gca,'xscale','linear','yscale','log','FontSize',8,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
end
end