function [FR_post_ref SP2]=plot_FR_bar(dat,xtick,cmap,y_lim,cat,t_mouse)

SP=[]; Arc_class=[];

for m=1:size(dat,2) %moue
    switch cat
        case 1 %cat 1: sigma dFF,
            for c1=1:size(dat{m}.Cal,1)
                for d1=1:3
                    if size(dat{m}.Cal{c1,d1})==1
                        dFF_sum(c1,d1)=NaN;
                    else
                        dFF_sum(c1,d1)=mean(dat{m}.Cal{c1,3+d1},'omitnan');
                    end
                end
            end
            SP=[SP;dFF_sum];
            ylab='Mean \DeltaF/F';
        case 2 % 2: sigma Peak,
            SP=[SP;cellfun(@extrack_sigma_peak,dat{m}.Peak)];
            ylab='\Sigma Peaks (\DeltaF/F)';
        case 3 % 3: Sigma transients
            for c1=1:size(dat{m}.Cal,1)
                for d1=1:3
                    Tr_area(c1,d1)=extract_area_dFF(dat{m}.Cal{c1,d1},dat{m}.Cal{c1,6+d1});
                end
            end
            SP=[SP;Tr_area];
            ylab='\Sigma Transients (\DeltaF/F)';
        case 4
            ylab='Ca^2^+ event rate (Hz)';
            %             [cellnumber length]=cellfun(@size,dat{m}.Cal(:,1:3));
            %             length=max(length,1);
            %             [SP2{m} rw]=cellfun(@size,dat{m}.Peak);
            %             SP2{m}(rw==1)=NaN;
            %             SP=[SP;SP2{m}./length*30];
            %             SP_tmp=SP2{m}./length*30;
            [length]=cellfun(@timelength,dat{m}.Cal(:,4:6));
            length=max(length,1);
            isnanPeriod=cellfun(@isnan,dat{m}.Cal(:,4:6),'UniformOutput',false);
            Peak=rmvnanPeak(dat{m}.Peak,isnanPeriod);
            [SP2{m} rw]=cellfun(@size,Peak);
            SP2{m}(rw==1)=NaN;
            SP=[SP;SP2{m}./length*30];
            SP_tmp=SP2{m}./length*30;
            Arc_tmp=dat{m}.Arc(:,2:end);
            if t_mouse
                for cls=1:4
                    FR{m,cls}=SP_tmp(find(Arc_tmp==cls));
                end
            end
        case 5 % Number of big transients
            ylab='Ca^2^+ event rate (Hz)';
            %SP{m}=cellfun(@decon_sum,dat{m}.S(:,1:3));
            SP2{m}=cellfun(@decon_sum,dat{m}.S(:,4:6));
            [cw rw]=cellfun(@size,dat{m}.S(:,1:3));
            SP2{m}(rw==1)=NaN;
            SP=[SP;SP2{m}];
    end
    Arc_class=[Arc_class;dat{m}.Arc(:,2:end)];
end

% for m=1:size(dat,2)
%     FR_cell{m,1}=SP_cell{m}(AC_cell{m}==1);
%     FR_cell{m,2}=SP_cell{m}(AC_cell{m}==2);
% end
SP=SP(:); Arc_class=Arc_class(:);
if t_mouse
    for mm=1:6
    FR_post_ref{mm,1}=SP2{mm}(find(dat{mm}.Arc_post_ref(:,2:end)==1));
    FR_post_ref{mm,2}=SP2{mm}(find(dat{mm}.Arc_post_ref(:,2:end)==2));
    end
else
    FR{1}=SP(find(Arc_class==1)); %No TXN Sum peak
    FR{2}=SP(find(Arc_class==2)); %TXN Sum peak
    FR{3}=SP(find(Arc_class==3));
    FR{4}=SP(find(Arc_class==4));
FR_post_ref{1}=[FR{1};FR{3}];
FR_post_ref{2}=[FR{2};FR{4}];
%     [p_anov anov_tbl stats]=anova1(cell2mat(FR'),[repmat(categorical({'~Pre ~Post'}),size(FR{1},1),1);...
%         repmat(categorical({'~Pre Post'}),size(FR{2},1),1);
%         repmat(categorical({'Pre ~Post'}),size(FR{3},1),1);...
%         repmat(categorical({'Pre Post'}),size(FR{4},1),1)]);
%     figure
%     [c,~,~,gnames] = multcompare(stats);
    [p a]=ranksum(FR_post_ref{1},FR_post_ref{2});
    M=[mean(FR_post_ref{1},'omitnan') mean(FR_post_ref{2},'omitnan')];
    S=[std(FR_post_ref{1},0,1,'omitnan')/sqrt(size(FR_post_ref{1},1)) std(FR_post_ref{2},0,1,'omitnan')/sqrt(size(FR_post_ref{2},1))];
    %S=[std(FR_post_ref{1},0,1,'omitnan') std(FR_post_ref{2},0,1,'omitnan')];
    figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
        'Color',[1 1 1],...
        'Renderer','painters','position',[100 100 300 300]);
    
    categ=[];
    for i=1:2
    categ=[categ; repmat(xtick{i},size(FR_post_ref{i},1),1)]; end
    %boxplot(cell2mat(FR_post_ref'),categ,'colors',cmap,'PlotStyle','compact')
    violin(FR_post_ref,'facecolor',cmap,'edgecolor',[],'bw',0.1,'mc',[],'medc',[],'plotlegend',0);
    hold all
    errorbar([1 2],M,S,'LineWidth',2,'linestyle','none','color','k','Capsize',10)
%     b=bar([1 2],M,'Barwidth',0.7,'LineWidth',2);
%     b.FaceColor='flat';
%     for i=1:2
%         b.CData(i,:)= cmap(i,:);
%         hold all
%         % plot(X,datum{i},'marker','.','markersize',2,'linestyle','none')
%     end
    set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:4],'XTickLabel',...
        xtick,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
    ylabel(ylab,'LineWidth',2,'FontSize',13,...
        'FontName','arial rounded mt bold')
    xlim([0.3 4.7])
    ylim([0 y_lim])
    p
    if p<0.05
        star='*';
        if p<0.01
            star='**';
            
            if p<0.001
                star='***';
                
            end
        end
    end
    if a
line([1 2],[y_lim*0.72+y_lim*0.05*i y_lim*0.72+y_lim*0.05*i],'color','k','linewidth',2);
    line([1 1],[y_lim*0.62+y_lim*0.05*i y_lim*0.72+y_lim*0.05*i],'color','k','linewidth',2);
    line([2 2],[y_lim*0.62+y_lim*0.05*i y_lim*0.72+y_lim*0.05*i],'color','k','linewidth',2);
        
        text(1.5,y_lim*(0.83),star,'HorizontalAlignment', 'center'...
            ,'LineWidth',2,'FontSize',13,'FontName','arial rounded mt bold')
    end
end


figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
        'Color',[1 1 1],'Renderer','painters','position',[100 100 300 300]);

bin=20;   cmap2=[0.1 0.1 0.1;1 0 0];
max_cdf=max(cellfun(@median_omitnan,FR_post_ref))*5;
for i=1:2
    
    [h{i} lag]=histcounts(FR_post_ref{i}(~isnan(FR_post_ref{i})),[0:max_cdf/bin:max_cdf],...
                   'Normalization','cdf');
plot(lag(1:end-1),h{i},'color',cmap2(i,:),'linewidth',2);  
hold all
end
xlabel('Ca^{2+} event rate (Hz)','FontSize',13,'FontName','arial rounded mt bold')
ylabel('CDF','FontSize',13,'FontName','arial rounded mt bold')
legend({'Arc^-','Arc^+',},'FontSize',13,'FontName','arial rounded mt bold')
set(gca,'linewidth',2,'FontSize',13,'FontName','arial rounded mt bold')
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
function M=median_omitnan(S)
M=median(S,'omitnan');
end