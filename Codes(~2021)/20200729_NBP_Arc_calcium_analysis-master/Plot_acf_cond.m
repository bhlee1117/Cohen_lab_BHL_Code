function acf=Plot_acf_cond(dat,groups,xtick,y_lim,x_lim,sw,bin,post_ref)
for d1=1:3
    for g=1:size(groups,2)
        acf{g,d1}=[];
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
        for l=find(list{m,g})'
            for d1=1:3
                
                if size(dat{m}.Cal{l,d1},2)>1 && ~isnan(dat{m}.Cal{l,d1}(1,1))
                    if isempty(bin) %No binning
                        datum=dat{m}.S{l,3+d1};
                        datum(dat{m}.run{d1}==0)=NaN;
                        [normacf lag]=autocorr(datum,'NumLags',1000);
                    else
                        [r c]=cellfun(@size,dat{m}.Cal(:,3+d1));  max_bin=ceil(max(c)*1/30/bin);
                        bin_tmp=zeros(1,max_bin);
                        if ~isempty(dat{m}.Peak{l,d1})
                            bin_tmp(1,ceil(dat{m}.Peak{l,d1}(:,3)/bin))=1;
                        end
                        nan_time=find(isnan(dat{m}.Cal{l,3+d1}))/30;
                        bin_tmp(1,ceil(nan_time/bin))=NaN; %out-of focus & grooming convert to NaN
                        [normacf lag]=autocorr(bin_tmp,'NumLags',round(1000/30/bin));
                    end
                    
                    acf{g,d1}=[acf{g,d1}; normacf];
                else
                    if isempty(bin)
                        acf{g,d1}=[acf{g,d1}; NaN(1,1001)];
                    else
                        acf{g,d1}=[acf{g,d1}; NaN(1,round(1000/30/bin)+1)];
                    end
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
            if isempty(bin)
                mseb([1:1:size(acf{g,d1},2)]*1/30,mean(acf{g,d1},1,'omitnan'),std(acf{g,d1},0,1,'omitnan')/sqrt(size(acf{g,d1},1)),lineProps,0.5)
            else
                mseb([1:1:size(acf{g,d1},2)]*bin,mean(acf{g,d1},1,'omitnan'),std(acf{g,d1},0,1,'omitnan')/sqrt(size(acf{g,d1},1)),lineProps,0.5)
            end
            hold all
        end
        ylabel('Auto-correlation')
        xlabel('Time (sec)')
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
            if isempty(bin)
                mseb([1:1:size(acf{g,d1},2)]*1/30,mean(acf{g,d1},1,'omitnan'),std(acf{g,d1},0,1,'omitnan')/sqrt(size(acf{g,d1},1)),lineProps,0.5)
            else
                mseb([0:1:size(acf{g,d1},2)-1]*bin,mean(acf{g,d1},1,'omitnan'),std(acf{g,d1},0,1,'omitnan')/sqrt(size(acf{g,d1},1)),lineProps,0.5)
            end
            hold all
        end
        ylabel('Auto-correlation')
        xlabel('Time (sec)')
        xlim([0 x_lim])
        ylim([0.001 y_lim])
        legend(xtick)
        set(gca,'yscale','log')
    end
end
end