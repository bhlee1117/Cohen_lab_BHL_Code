function [FR]=plot_dFFarea_Cond(dat,groups,xtick,ylab,y_lim,sw)

for m=1:size(dat,2) %mouse
    for c1=1:size(dat{m}.Cal,1)
        for d1=1:3
       if isnan(dat{m}.Cal{c1,d1})
SP{m}(c1,d1)=NaN;
       else
%SP{m}(c1,d1)=extract_area_dFF(dat{m}.Cal{c1,d1},dat{m}.Cal{c1,6+d1});
SP{m}(c1,d1)=sum(dat{m}.Cal{c1,d1});
       end
        end
    end
Arc_class{m}=[dat{m}.Arc(:,2:end)];
for g=1:size(groups,2)
for cd=1:size(groups{g},2)/2
    list{m,g}(:,cd)=Arc_class{m}(:,groups{g}(1,2*cd-1))==groups{g}(1,2*cd);
end
list{m,g}=sum(list{m,g},2)==size(groups{g},2)/2; % cells that match the conditions
FR{m,g}=[SP{m}(find(list{m,g}),:)];
end
end
p_list=[1 2;1 3;2 3];
if size(groups,2)==2
p_list=p_list(1,:);
end


if sw
for g=1:size(groups,2)
    figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 300 300]);
d=cell2mat(FR(:,g));
for pp=1:size(p_list,1)
[a p_value(pp)]=ttest(d(:,p_list(pp,1)),d(:,p_list(pp,2))); end
p_value
M=mean(d,1,'omitnan');
n=sum(~isnan(d),1);
S=std(d,0,1,'omitnan')./sqrt(n);
G=[1 2 3];
errorbar(G,M,S,...
        'LineWidth',2,'linestyle','-','color','k','Capsize',10)
    hold all
% b=bar(G,M,'Barwidth',0.7,'LineWidth',2);
% b.FaceColor='flat';
% for i=1:2
% b.CData(i,:)= cmap(i,:);
% hold all
% % plot(X,datum{i},'marker','.','markersize',2,'linestyle','none')
% end
set(gca,'FontSize',8,'LineWidth',1,'XTick',G,'XTickLabel',...
 xtick,'FontName','arial rounded mt bold','FontSize',13,'LineWidth',2);
ylabel(ylab,'LineWidth',2,'FontSize',13,...
    'FontName','arial rounded mt bold')
xlim([0.3 3.7])
ylim([0 y_lim])
for i=1:3
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
   
   for day=1:3
       for g=1:size(groups,2)
tmp=cell2mat(FR(:,g));
d{day,g}=tmp(:,day);
       end
figure1 = figure('InvertHardcopy','off','PaperUnits','centimeters',...
    'Color',[1 1 1],...
    'Renderer','painters','position',[100 100 300 300]);
for pp=1:size(p_list,1)
[a p_value(pp)]=kstest2(d{day,p_list(pp,1)},d{day,p_list(pp,2)}); 
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
for i=1:size(p_list,1)
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
end
end
