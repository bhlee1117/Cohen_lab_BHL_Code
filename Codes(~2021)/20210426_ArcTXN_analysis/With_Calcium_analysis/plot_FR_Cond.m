function FR=plot_FR_Cond(dat,groups,xtick,y_lim,sw,cat,post_ref,t_mouse)

for m=1:size(dat,2) %mouse
    switch cat
        case 1 %Mean dFF
            ylab='Mean \DeltaF/F';
            for c1=1:size(dat{m}.Cal,1)
                for d1=1:3
                    if size(dat{m}.Cal{c1,d1},2)<2
                        SP{m}(c1,d1)=NaN;
                    else
                        tr=dat{m}.Cal{c1,3+d1};
                        tr(dat{m}.run{d1}==0)=NaN;
                        SP{m}(c1,d1)=mean(tr,'omitnan');
                    end
                end
            end
        case 2 % Sum peak
            ylab='\Sigma Peaks';
            clear P
            %             for c1=1:size(dat{m}.Cal,1)
            %                 for d1=1:3
            %                     if size(dat{m}.Peak{c1,d1},2)<2
            %                     P(c1,d1)=NaN;
            %                     else if
            %                     dat{m}.Peak{c1,d1}(find(dat{m}.run{d1}(round(dat{m}.Peak{c1,d1}(:,3)*30))==0),:)=[];
            %                     P(c1,d1)=sum(dat{m}.Peak{c1,d1}(:,2))/sum(dat{m}.run{d1}==1 & ~isnan(dat{m}.Cal{c1,3+d1}))*30;
            %                 end
            %             end
            SP{m}=cellfun(@extrack_sigma_peak,dat{m}.Peak);
        case 3 % Sum Transient
            ylab='\Sigma Transient (\DeltaF/F)';
            for c1=1:size(dat{m}.Cal,1)
                for d1=1:3
                    if size(dat{m}.Cal{c1,d1},2)<2
                        SP{m}(c1,d1)=NaN;
                    else
                        tr_bin=zeros(size(dat{m}.run{d1},1),1);
                        dat{m}.Cal{c1,6+d1}(find(dat{m}.Cal{c1,6+d1}(:,1)>size(dat{m}.Cal{c1,3+d1},2)),:)=[];
                        for tr=1:size(dat{m}.Cal{c1,6+d1},1)
                            if dat{m}.Cal{c1,6+d1}(tr,2)>size(dat{m}.run{d1},1)
                                dat{m}.Cal{c1,6+d1}(tr,2)=size(dat{m}.run{d1},1);
                            end
                            if sum(dat{m}.Cal{c1,6+d1}(tr,:)==0)==0 %
                                tr_bin(dat{m}.Cal{c1,6+d1}(tr,1):dat{m}.Cal{c1,6+d1}(tr,2),1)=1;
                            end
                        end
                        tr_bin(dat{m}.run{d1}==0)=NaN;  tr_bin(isnan(dat{m}.Cal{c1,3+d1}))=NaN;
                        SP{m}(c1,d1)=sum(dat{m}.Cal{c1,d1}(1,find(tr_bin==1)))/sum(~isnan(tr_bin))*30;
                    end
                end
            end
        case 4 % Number of peaks
            ylab='Ca^2^+ event rate (Hz)';
            [cellnumber length]=cellfun(@size,dat{m}.Cal(:,1:3));
            length=max(length,1);
            [SP{m} rw]=cellfun(@size,dat{m}.Peak);
            SP{m}(rw==1)=NaN;
            SP{m}=SP{m}./length*30;
        case 5 % Number of big transients
             ylab='Ca^2^+ event rate (Hz)';
             %SP{m}=cellfun(@decon_sum,dat{m}.S(:,1:3));
             SP{m}=cellfun(@decon_sum,dat{m}.S(:,4:6));
             [cw rw]=cellfun(@size,dat{m}.Peak);
             SP{m}(rw==1)=NaN;
        case 6
            ylab='Maximum firing rate';
            for c1=1:size(dat{m}.S,1)
                for d1=1:3
            tmp{c1,d1}=movsum(dat{m}.S{c1,3+d1},30,'omitnan');
                end
            end
            %SP{m}=cellfun(@max,dat{m}.S(:,4:6));
            [cw rw]=cellfun(@size,dat{m}.Peak);
            SP2{m}(rw==1)=NaN;
            SP{m}=cellfun(@max,tmp);
            
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
        FR{m,g}=[SP{m}(find(list{m,g}),:)];
        % Matrix FR : cell array consist of n X 3 matrix, each column represents
        % days.
        % Quantity of interest is arranged in mouse and groups.
        
    end
end
p_list=[1 2;2 3;1 3;1 4;2 4;3 4];

%%%%% Data arrange done %%%%%
%%%%% Significance check & plotting start %%%%%

if sw
    p_list=p_list(1:3,:);
    for g=1:size(groups,2)
        figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
            'Color',[1 1 1],...
            'Renderer','painters','position',[100 100 300 300]);
        if t_mouse
            for mm=1:size(FR,1)
                d(mm,:)=mean(FR{mm,g},'omitnan');
            end
        else
            d=cell2mat(FR(:,g));
        end
        for pp=1:size(p_list,1)
            %[a p_value(pp)]=ttest(d(:,p_list(pp,1)),d(:,p_list(pp,2)));
            [p_value(pp) a]=ttest(d(:,p_list(pp,1)),d(:,p_list(pp,2)));
        end
        M=mean(d,1,'omitnan');
        n=sum(~isnan(d),1);
        S=std(d,0,1,'omitnan')./sqrt(n);
        G=[0.7 1 1.3];
        plot_errorbar2(G,d,'ttest',y_lim,ylab,xtick,0)
%         for i=1:size(d,1)
%             plot(G,d(i,:),'color',[0.7 0.7 0.7],'marker','.')
%             hold all
%         end
%         errorbar(G,M,S,...
%             'LineWidth',2,'linestyle','-','color','k','Capsize',10)
%         hold all
%         set(gca,'FontSize',8,'LineWidth',1,'XTick',G,'XTickLabel',...
%             xtick,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
%         ylabel(ylab,'LineWidth',2,'FontSize',13,...
%             'FontName','arial rounded mt bold')
%         xlim([0.3 3.7])
%         ylim([0 y_lim])
    end
else  % Day -> groups
    if size(groups,2)==2
        p_list=p_list(1,:); end
    if size(groups,2)==3
        p_list=p_list(1:3,:); end
    if size(groups,2)==4
        p_list=p_list(1:6,:); end
    
    figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
        'Color',[1 1 1],...
        'Renderer','painters','position',[100 100 250+size(groups,2)*50 250]);
    if size(groups,2)==2
        G=[0.8 1.2 1.8 2.2 2.8 3.2];
    else if size(groups,2)==3
        G=[0.7 1 1.3 1.7 2 2.3 2.7 3 3.3];
        else
        G=[0.7 0.9 1.1 1.3 1.7 1.9 2.1 2.3 2.7 2.9 3.1 3.3];
        end
    end
    for day=1:3
        if t_mouse
            for g=1:size(groups,2)
                for kk=1:size(FR,1)
                    tmp(kk,:)=mean(FR{kk,g},1,'omitnan'); end
                d{day,g}=tmp(:,day);
            end
        else
            for g=1:size(groups,2)
                tmp=cell2mat(FR(:,g));
                d{day,g}=tmp(:,day);
            end
        end
        M=[]; S=[];
        cmap=[0.1 0.1 0.1;distinguishable_colors(3)]; gn=size(groups,2);
        for g=1:size(groups,2)
            M=[M mean(d{day,g},1,'omitnan')];
            S=[S  std(d{day,g},0,1,'omitnan')./sqrt(sum(~isnan(d{day,g})))]; end
        M
        b=bar(G(gn*(day-1)+1:gn*day),M,'Barwidth',0.7,'LineWidth',2);
        hold all
        b.FaceColor='flat';
        for j=1:size(groups,2)
            b.CData(j,:)= cmap(j,:);
        end
        
        if t_mouse
            for pp=1:size(p_list,1)
                [a p_value(pp,1)]=ttest(d{day,p_list(pp,1)},d{day,p_list(pp,2)});
            end
            for mice=1:size(d{day,g},1)
                tmp2=cell2mat(d(day,:));
                %plot([1:size(groups,2)]+randn(1,1)*0,tmp2(mice,:),'marker','o','color',[0.4 0.4 0.4],'linewidth',1)
                plot(G([gn*(day-1)+1:gn*day]),tmp2(mice,:),'marker','o','color',[0.7 0.7 0.7],'linewidth',1)
                hold all
            end
        else
            for pp=1:size(p_list,1)
                [a p_value(pp,1)]=kstest2(d{day,p_list(pp,1)},d{day,p_list(pp,2)});
                [a p_value(pp,2)]=ttest2(d{day,p_list(pp,1)},d{day,p_list(pp,2)});
                [p_value(pp,3) a]=ranksum(d{day,p_list(pp,1)},d{day,p_list(pp,2)});
            end
            for g=1:size(groups,2)
                plot(G(gn*(day-1)+g)+randn(size(d{day,g},1),1)*0.01,d{day,g},'marker','.','color',[0.7 0.7 0.7],'linestyle','none')
                hold all
            end
        end
        p_value;
        errorbar(G(gn*(day-1)+1:gn*day),M,S,...
            'LineWidth',2,'linestyle','-','color','k','Capsize',10)
        
        set(gca,'FontSize',8,'LineWidth',1,'XTick',G,'XTickLabel',...
            xtick,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
        ylabel(ylab,'LineWidth',2,'FontSize',13,...
            'FontName','arial rounded mt bold')
        %xlim([0.3 size(groups,2)+0.7])
        xlim([0.5 G(end)+0.3])
        ylim([0 y_lim])
        if ~t_mouse
            ref=1;
        else
            ref=1; end
        for i=1:size(p_list,1)
            if p_value(i,ref)<0.05
                star='*';
                if p_value(i,ref)<0.01
                    star='**';
                    if p_value(i,ref)<0.001
                        star='***';
                    end
                end
                line(day-1+G(p_list(i,:)),[y_lim*0.77+y_lim*0.05*i y_lim*0.77+y_lim*0.05*i],'color','k','linewidth',2);
                line(day-1+G([p_list(i,1) p_list(i,1)]),[y_lim*0.70+y_lim*0.05*i y_lim*0.77+y_lim*0.05*i],'color','k','linewidth',2);
                line(day-1+G([p_list(i,2) p_list(i,2)]),[y_lim*0.70+y_lim*0.05*i y_lim*0.77+y_lim*0.05*i],'color','k','linewidth',2);
                text(day-1+mean(G(p_list(i,:))),y_lim*0.79+y_lim*0.05*i,star,'FontSize',13,'FontName','arial rounded mt bold',...
                    'HorizontalAlignment', 'center')
            end
        end
        xtickangle(45)
    end
end

end
function length=timelength(tr)
length=sum(~isnan(tr));
end
function newpeak=rmvnanPeak(Peak,isnanPeriod)
newpeak=Peak;
for c=1:size(Peak,1)
    for d1=1:size(Peak,2)
        if size(Peak{c,d1},2)>1
            newpeak{c,d1}(find(isnanPeriod{c,d1}(round(Peak{c,d1}(:,3)*30))),:)=[];
        end
    end
end
end

