function [HR]=plot_histo_bin_trace(dat,groups,tick,sw,bin,post_ref,t_mouse)

for m=1:size(dat,2) %mouse
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
        BR{m,g}=[binned_tr{m}(find(list{m,g}),:)];
    end
end
p_list=[1 2;1 3;2 3;1 4;2 4;3 4];
% Now plot
if t_mouse
    
    if sw
        cmap=distinguishable_colors(3);
        for g=1:size(groups,2)
            HR{g,d1}=[];
            figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
                'Color',[1 1 1],'Renderer','painters','position',[100 100 300 200]);
            
            for d1=1:3
                for m=1:size(dat,2)
                    h_tmp=[];
                    for c=1:size(BR{m,g},1)
                        %fig=figure; set(fig,'visible','off');
                        %h=histcounts(BR{m,g}{c,d1}(find(~isnan(BR{m,g}{c,d1}))),[0:1:30],'Normalization','Probability');
                        h=histcounts(BR{m,g}{c,d1}(find((BR{m,g}{c,d1})>0)),[0:1:30],'Normalization','Probability');
                        h_tmp=[h_tmp; h];
                        %close(fig)
                    end
                    
                    HR{g,d1}=[HR{g,d1}; mean(h_tmp,1,'omitnan')];
                end
                lineProps.col{1}=cmap(d1,:);
                mseb([0:size(HR{g,d1},2)-1]*bin,mean(HR{g,d1},1,'omitnan')...
                    ,std(HR{g,d1},0,1,'omitnan')/sqrt(size(HR{g,d1},1)),lineProps,1);
                hold all
            end
            legend(tick)
            xlabel(ylab,'LineWidth',2,'FontSize',13,...
                'FontName','arial rounded mt bold')
            ylabel('Probability','LineWidth',2,'FontSize',13,...
                'FontName','arial rounded mt bold')
            set(gca,'xscale','linear','yscale','log','FontSize',8,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
        end
        
    else
        cmap=[0.2 0.2 0.2;1 0 0;1 0 1;0 0 1];
        for d1=1:3
            HR{g,d1}=[];
            figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
                'Color',[1 1 1],'Renderer','painters','position',[100 100 300 200]);
            
            for g=1:size(groups,2)
                for m=1:size(dat,2)
                    h_tmp=[];
                    for c=1:size(BR{m,g},1)
                        %fig=figure; set(fig,'visible','off');
                        %h=histcounts(BR{m,g}{c,d1}(find(~isnan(BR{m,g}{c,d1}))),[0:1:30],'Normalization','Probability');
                        h=histcounts(BR{m,g}{c,d1}(find((BR{m,g}{c,d1})>0)),[0:1:30],'Normalization','Probability');
                        h_tmp=[h_tmp; h];
                        %close(fig)
                    end
                    HR{g,d1}=[HR{g,d1}; mean(h_tmp,1,'omitnan')];
                end
                lineProps.col{1}=cmap(g,:);
                mseb([0:size(HR{g,d1},2)-1]*bin,mean(HR{g,d1},1,'omitnan')...
                    ,std(HR{g,d1},0,1,'omitnan')/sqrt(size(HR{g,d1},1)),lineProps,1);
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
    
else
    
    if sw
        cmap=distinguishable_colors(3);
        for g=1:size(groups,2)
            HR{g,d1}=[];
            figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
                'Color',[1 1 1],'Renderer','painters','position',[100 100 300 200]);
            
            for d1=1:3
                for m=1:size(dat,2)
                    for c=1:size(BR{m,g},1)
                        %fig=figure; set(fig,'visible','off');
                        %h=histcounts(BR{m,g}{c,d1}(find(~isnan(BR{m,g}{c,d1}))),[0:1:30],'Normalization','Probability');
                        h=histcounts(BR{m,g}{c,d1}(find((BR{m,g}{c,d1})>0)),[0:1:30],'Normalization','Probability');
                        HR{g,d1}=[HR{g,d1}; h];
                        %close(fig)
                    end
                end
                lineProps.col{1}=cmap(d1,:);
                mseb([0:size(HR{g,d1},2)-1]*bin,mean(HR{g,d1},1,'omitnan')...
                    ,std(HR{g,d1},0,1,'omitnan')/sqrt(size(HR{g,d1},1)),lineProps,1);
                hold all
            end
            legend(tick)
            xlabel(ylab,'LineWidth',2,'FontSize',13,...
                'FontName','arial rounded mt bold')
            ylabel('Probability','LineWidth',2,'FontSize',13,...
                'FontName','arial rounded mt bold')
            set(gca,'xscale','linear','yscale','log','FontSize',8,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
        end
        
    else
        cmap=[0.2 0.2 0.2;1 0 0];
        for d1=1:3
            HR{g,d1}=[];
            figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
                'Color',[1 1 1],'Renderer','painters','position',[100 100 300 200]);
            
            for g=1:size(groups,2)
                for m=1:size(dat,2)
                    for c=1:size(BR{m,g},1)
                        %fig=figure; set(fig,'visible','off');
                        %h=histcounts(BR{m,g}{c,d1}(find(~isnan(BR{m,g}{c,d1}))),[0:1:30],'Normalization','Probability');
                        h=histcounts(BR{m,g}{c,d1}(find((BR{m,g}{c,d1})>0)),[0:1:30],'Normalization','Probability');
                        HR{g,d1}=[HR{g,d1}; h];
                        %close(fig)
                    end
                end
                lineProps.col{1}=cmap(g,:);
                mseb([0:size(HR{g,d1},2)-1]*bin,mean(HR{g,d1},1,'omitnan')...
                    ,std(HR{g,d1},0,1,'omitnan')/sqrt(size(HR{g,d1},1)),lineProps,1);
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
end
end
function tr=peak_to_trace(peak,cal_tr,bin,maxbin)
tr=zeros(1,size(cal_tr,2));
%tr=zeros(1,maxbin);
if size(peak,2)>1
    %         [C,ia,ic] = unique(ceil(peak(:,3)/bin));
    %         a_counts = accumarray(ic,1);
    %     tr(1,C)=a_counts;
    tr(1,round(peak(:,3)*30))=1;
end
tr(1,find(isnan(cal_tr)))=NaN;
%tr(1,ceil(find(isnan(cal_tr))/30))=NaN;
tr=movsum(tr,bin*30);
end