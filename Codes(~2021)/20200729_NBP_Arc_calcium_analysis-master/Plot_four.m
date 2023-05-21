function F=Plot_four(dat,xtick,y_lim,x_lim,post_ref)
if post_ref
    cd_max=2;
else
    cd_max=4;
end
for cd=1:cd_max
    F{cd,1}=[];  end

for m=1:size(dat,2)
    if post_ref
    Arc_class{m}=[dat{m}.Arc_post_ref(:,2:end)];
    else
    Arc_class{m}=[dat{m}.Arc(:,2:end)];
    end
    datum=dat{m}.FT(:,4:6);
    datum=datum(:);
    
    
    for cd=1:cd_max
        list{m,cd}=Arc_class{m}==cd;
        list{m,cd}= list{m,cd}(:);
        for c=find(list{m,cd})'
            if size(datum{c},2)==1
                F{cd,1}=[F{cd,1}; NaN(1,6001)];
            else
                [normF]=datum{c};
                F{cd,1}=[F{cd,1}; normF];
            end
        end
    end
end

%% plot

cmap=distinguishable_colors(4);

figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 400 250]);

for g=1:cd_max
    lineProps.col{1}=cmap(g,:);
    mseb([0:15/6000:15],mean(F{g,1},1,'omitnan'),std(F{g,1},0,1,'omitnan')./sqrt(sum(~isnan(F{g,1}),1)),lineProps,0.5);
hold all
end


ylabel('Power')
xlabel('Frequency (Hz)')
xlim([0 x_lim])
ylim([0.001 y_lim])
legend(xtick)
set(gca,'yscale','log')
end

