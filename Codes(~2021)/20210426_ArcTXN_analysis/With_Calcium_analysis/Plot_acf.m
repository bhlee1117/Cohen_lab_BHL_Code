function acf=Plot_acf(dat,xtick,y_lim,x_lim)
for cd=1:4
    acf{cd,1}=[];  end

for m=1:size(dat,2)
    Arc_class{m}=[dat{m}.Arc(:,2:end)];
    
    datum=dat{m}.S(:,4:6);
    datum=datum(:);
    
    
    for cd=1:4
        list{m,cd}=Arc_class{m}==cd;
        list{m,cd}= list{m,cd}(:);
        for c=find(list{m,cd})'
            if size(datum{c},2)==1
                acf{cd,1}=[acf{cd,1}; NaN(1,1001)];
            else
                [normacf lag]=autocorr(datum{c},'NumLags',1000);
                acf{cd,1}=[acf{cd,1}; normacf];
            end
        end
    end
end

%% plot

cmap=distinguishable_colors(4);

figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 400 250]);

for g=1:4
    lineProps.col{1}=cmap(g,:);
    mseb([1:1:size(acf{g,1},2)]*1/30,mean(acf{g,1},1,'omitnan'),std(acf{g,1},0,1,'omitnan')./sqrt(sum(~isnan(acf{g,1}),1)),lineProps,0.5)
hold all
end


ylabel('Auto-correlation')
xlabel('Time (sec)')
xlim([0 x_lim])
ylim([0.001 y_lim])
legend(xtick)
set(gca,'yscale','log')
end

