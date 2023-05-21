%% Load mouse data

clear
[fnm pth]=uigetfile('.mat','Select the analyzed mouse data','Multiselect','on');

%%
clear CFC HC
if ischar(fnm)
    fnm={fnm};
end
lg=[];
h=1; r=1;
for i=1:length(fnm)
    load([pth fnm{1,i}])
    sp=split(fnm{1,i},'_');
    switch sp{3,1}(1)
        case 'C'
    
    CFC{h}=RNAscope.result(find(sum(RNAscope.result(:,3:5)==3,2)==0),:);
    h=h+1;
        case 'H'
          HC{r}=RNAscope.result(find(sum(RNAscope.result(:,3:5)==3,2)==0),:);
          r=r+1;
    end
end

%% 'Fraction of TXN Cell (%)'
clear SEM
    figure(1)
   p=[];

for i=1:size(CFC,2)
    for c=1:3
    fracTXN{1}(i,c)=sum(CFC{1,i}(:,2+c)==2)/sum(CFC{1,i}(:,2+c)~=3);
    end
end
SEM(1,:)=std(fracTXN{1},0,1,'omitnan')/sqrt(size(CFC,2));

for i=1:size(HC,2)
    for c=1:3
        fracTXN{2}(i,c)=sum(HC{1,i}(:,2+c)==2)/sum(HC{1,i}(:,2+c)~=3);
    end
end
SEM(2,:)=std(fracTXN{2},0,1,'omitnan')/sqrt(size(HC,2));

cc=[0.5 0.9 0.5;0.9 0.5 0.5;0.5 0.5 0.9];
b=bar([1 2],[mean(fracTXN{1},1) ; mean(fracTXN{2},1)]*100,'BarWidth',0.8);%'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
for j=1:3
b(j).FaceColor=cc(mod(j,3)+3*floor(j/3),1);
%b(j).CData= cc(mod(j,3)+3*floor(j/3),1);
end
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1 2],'XTickLabel',...
    {'CFC','HC'});
   ylabel('Fraction of TXN Cell (%)')
   xlim([0.5 2.5])
   hold all
   errorbar([0.77 1 1.23 1.77 2 2.23],[mean(fracTXN{1},1) mean(fracTXN{2},1)]*100,[SEM(1,:)*100 SEM(2,:)*100],'k','linewidth',4,'linestyle','none')
%xtickangle(45);
%%
color=[0 1 0;1 0 0;0 1 1];
gr=0.23;
for i=1:2 %% plot by line
    for j=1:size(fracTXN{i},1)
        plot([i-gr i i+gr],fracTXN{i}(j,1:end)*100,'.','markersize',15,'linestyle','-')
        hold all
    end
end
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1-gr 1 1+gr 2-gr 2 2+gr],'XTickLabel',...
    {'Arc CFC','C-fos CFC','Egr1 CFC','Arc HC','C-fos HC','Egr1 HC'});
   ylabel('Fraction of TXN Cell (%)')
   xtickangle(45);
%%
figure(2)

clear overlap over_name overlapTXN SEM data
Window_size=200;
exhib=[3 4 5];
for i=[1]
%HC CFC Ret1
  handles.fig = figure(1);
set(gcf, 'Position',  [50, 30, Window_size*(size(CFC,2)+2), 300])
lg=[{'CFC-2-1'},{'CFC-2-2'},{'CFC-2-3'},{'CFC-3-1'},{'CFC-3-2'},{'CFC-3-3'}];
for j=1:size(CFC,2)
  
handles.axes1 = axes('Units','pixels','Position',[Window_size*(j-1)+50 50 Window_size Window_size]);   
title(lg{i,j})
tmp=venn_data(CFC{1,j}(:,exhib))/sum(sum(CFC{1,j}==2,2)>0);%size(CFC{1,j},1);
data{i}(j,:)=tmp;
[a b]=venn(data{i}(j,:));
for t=1:size(data{i},2)
    text(b.ZoneCentroid(t,1)-b.ZoneArea(1,t)*0.1,b.ZoneCentroid(t,2),[num2str(data{i}(j,t)*100,2),' % '])
end
axis off equal tight
end
handles.axes1 = axes('Units','pixels','Position',[Window_size*(j)+100 50 Window_size Window_size]);   
title('Averaged')
[a b]=venn(mean(data{i},1));
ave_data=mean(data{i});
for t=1:size(data{i},2)
    text(b.ZoneCentroid(t,1)-b.ZoneArea(1,t)*0.1,b.ZoneCentroid(t,2),[num2str(ave_data(1,t)*100,2),' % '])
end
axis off equal tight
%Ret1, Ret2, Ret3
%CFC on, 123
%CFC off, 123
end
%%
exhib=[3:5];
for i=1:size(CFC,2)
    overlap{1}(i,:)=venn_data(CFC{1,i}(:,exhib));
end

for i=1:size(HC,2)
    overlap{2}(i,:)=venn_data(HC{1,i}(:,exhib));
end

for i=1:2
    for j=1:size(overlap{i},1)
        cond_prob{i}(j,:)=[(overlap{i}(j,4)+overlap{i}(j,7))/(overlap{i}(j,1)+overlap{i}(j,4)+overlap{i}(j,7)+overlap{i}(j,5)) overlap{i}(j,5)/overlap{i}(j,1) overlap{i}(j,6)/overlap{i}(j,2)...
                   overlap{i}(j,7)/overlap{i}(j,1) overlap{i}(j,7)/overlap{i}(j,2) overlap{i}(j,7)/overlap{i}(j,3)]; %c-fos&Arc|Arc, Arc&Egr-1|Arc, Egr-1&c-fos|c-fos, All|Arc, All|c-fos, All|Egr-1
        if i==1
            plot([1:6],cond_prob{i}(j,:),'.','markersize',15,'color',[0.8 0.7 0.7],'linestyle','-')
        else
            plot([1:6],cond_prob{i}(j,:),'.','markersize',15,'color',[0.7 0.8 0.7],'linestyle','-')
        end
        hold all
   set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:6],'XTickLabel',...
    {'C-fos & Arc | Arc', 'Arc & Egr-1 | Arc', 'Egr-1 & C-fos | C-fos', 'All | Arc', 'All | C-fos', 'All | Egr-1'});
   ylabel('Conditional probability (%)')
   xtickangle(45);
   xlim([0.5 6.5])
   
    end
     if i==1 
         plot([1:6],mean(cond_prob{i},1),'.','markersize',15,'color',[1 0 0],'linestyle','-')
     else
          plot([1:6],mean(cond_prob{i},1),'.','markersize',15,'color',[0 1 0],'linestyle','-')
     end
end