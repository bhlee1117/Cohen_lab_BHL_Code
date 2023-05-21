function FR=plot_FR_Cond_day(dat,groups,doi,xtick,y_lim,sw,cat)

for m=1:size(dat,2) %mouse
    switch cat
        case 1 %Mean dFF
            ylab='Mean \DeltaF/F';
            for c1=1:size(dat{m}.Cal,1)
                for d1=1:3
                    if size(dat{m}.Cal{c1,d1},2)<2
                        SP{m}(c1,d1)=NaN;
                    else
                        SP{m}(c1,d1)=mean(dat{m}.Cal{c1,3+d1},'omitnan');
                    end
                end
            end
        case 2 % Sum peak
            ylab='\Sigma Peaks';
            SP{m}=cellfun(@extrack_sigma_peak,dat{m}.Peak);
        case 3 % Sum Transient
            ylab='\Sigma Transient (\DeltaF/F)';
            for c1=1:size(dat{m}.Cal,1)
                for d1=1:3
                    if size(dat{m}.Cal{c1,d1},2)<2
                        SP{m}(c1,d1)=NaN;
                    else
                        tr_area=0;
                        for tr=1:size(dat{m}.Cal{c1,6+d1},1)
                            if sum(dat{m}.Cal{c1,6+d1}(tr,:)==0)==0
                            tr_area=tr_area+sum(dat{m}.Cal{c1,3+d1}(dat{m}.Cal{c1,6+d1}(tr,1):dat{m}.Cal{c1,6+d1}(tr,2)));
                            end
                        end
                        SP{m}(c1,d1)=tr_area;
                    end
                end
            end
        case 4 % Number of peaks
            ylab='Number of Ca^2^+ events';
            [SP{m} rw]=cellfun(@size,dat{m}.Peak);
            SP{m}(rw==1)=NaN;
        case 5 % Number of big transients
    end
    
    Arc_class{m}=[dat{m}.Arc(:,2:end)];
    for g=1:size(groups,2)
        for cd=1:size(groups{g},2)/2
            list{m,g}(:,cd)=Arc_class{m}(:,groups{g}(1,2*cd-1))==groups{g}(1,2*cd);
        end
        list{m,g}=sum(list{m,g},2)==size(groups{g},2)/2; % cells that match the conditions
        FR{m,g}=[SP{m}(find(list{m,g}),:)];
        % Matrix FR : cell array consist of n X 3 matrix, each column represents
        % days.
        % Quantity of interest is arranged in mouse and groups.
        
    end
end
p_list=[1 2;1 3;2 3;1 4;2 4;3 4];

%%%%% Data arrange done %%%%%
%%%%% Significance check & plotting start %%%%%

if sw
    p_list=p_list(1:3,:);
    for g=1:size(groups,2)
        figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
            'Color',[1 1 1],...
            'Renderer','painters','position',[100 100 300 300]);
        d=cell2mat(FR(:,g));
        for pp=1:size(p_list,1)
            [a p_value(pp)]=ttest(d(:,p_list(pp,1)),d(:,p_list(pp,2))); end
        M=mean(d(:,doi),1,'omitnan');
        n=sum(~isnan(d(:,doi)),1);
        S=std(d(:,doi),0,1,'omitnan')./sqrt(n);
        G=[1 2 3];
        for i=1:size(d,1)
            plot(G,d(i,doi),'color',[0.4 0.4 0.4],'marker','.')
            hold all
        end
        errorbar(G,M,S,...
            'LineWidth',2,'linestyle','-','color','k','Capsize',10)
        hold all
        set(gca,'FontSize',8,'LineWidth',1,'XTick',G,'XTickLabel',...
            xtick,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
        ylabel(ylab,'LineWidth',2,'FontSize',13,...
            'FontName','arial rounded mt bold')
        xlim([0.3 size(doi,2)+0.7])
        ylim([0 y_lim])
        
        for i=1:size(p_value,2)
            if p_value(i)<0.05
                star='*';
                if p_value(i)<0.01
                    star='**';
                    if p_value(i)<0.001
                        star='***';
                    end
                end
                
                line(p_list(i,:),[y_lim*0.8+y_lim*0.05*i y_lim*0.8+y_lim*0.05*i],'color','k','linewidth',2);
                line([p_list(i,1) p_list(i,1)],[y_lim*0.7+y_lim*0.05*i y_lim*0.8+y_lim*0.05*i],'color','k','linewidth',2);
                line([p_list(i,2) p_list(i,2)],[y_lim*0.7+y_lim*0.05*i y_lim*0.8+y_lim*0.05*i],'color','k','linewidth',2);
                text(mean(p_list(i,:)),y_lim*0.82+y_lim*0.05*i,star,'FontSize',13,'FontName','arial rounded mt bold',...
                    'HorizontalAlignment', 'center')
            end
        end
    end
else  % Day -> groups
    if size(groups,2)==2
    p_list=p_list(1,:); end
    if size(groups,2)==3
    p_list=p_list(1:3,:); end
    if size(groups,2)==4
    p_list=p_list(1:6,:); end
    for day=doi
        for g=1:size(groups,2)
            tmp=cell2mat(FR(:,g));
            d{day,g}=tmp(:,day);
        end
        figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
            'Color',[1 1 1],...
            'Renderer','painters','position',[100 100 300 300]);
        for pp=1:size(p_list,1)
            [a p_value(pp,1)]=kstest2(d{day,p_list(pp,1)},d{day,p_list(pp,2)});
            [a p_value(pp,2)]=ttest2(d{day,p_list(pp,1)},d{day,p_list(pp,2)});
            [p_value(pp,3) a]=ranksum(d{day,p_list(pp,1)},d{day,p_list(pp,2)});
        end
        p_value
        M=[]; S=[];
        for g=1:size(groups,2)
            M=[M mean(d{day,g},1,'omitnan')];
            S=[S  std(d{day,g},0,1,'omitnan')./sqrt(sum(~isnan(d{day,g})))]; end
        G=[1:size(groups,2)];
        
        errorbar(G,M,S,...
            'LineWidth',2,'linestyle','-','color','k','Capsize',10)
        hold all
        set(gca,'FontSize',8,'LineWidth',1,'XTick',G,'XTickLabel',...
            xtick,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
        ylabel(ylab,'LineWidth',2,'FontSize',13,...
            'FontName','arial rounded mt bold')
        xlim([0.3 size(groups,2)+0.7])
        ylim([0 y_lim])
        ref=1;
        for i=1:size(p_list,1)
            if p_value(i,ref)<0.05
                star='*';
                if p_value(i,ref)<0.01
                    star='**';
                    if p_value(i,ref)<0.001
                        star='***';
                    end
                end
                
                line(p_list(i,:),[y_lim*0.72+y_lim*0.04*i y_lim*0.72+y_lim*0.04*i],'color','k','linewidth',2);
                line([p_list(i,1) p_list(i,1)],[y_lim*0.62+y_lim*0.04*i y_lim*0.72+y_lim*0.04*i],'color','k','linewidth',2);
                line([p_list(i,2) p_list(i,2)],[y_lim*0.62+y_lim*0.04*i y_lim*0.72+y_lim*0.04*i],'color','k','linewidth',2);
                text(mean(p_list(i,:)),y_lim*0.74+y_lim*0.04*i,star,'FontSize',13,'FontName','arial rounded mt bold',...
                    'HorizontalAlignment', 'center')
            end
        end
    end
end
end
