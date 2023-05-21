%% Load mouse data

clear
[fnm pth]=uigetfile('.mat','Select the analyzed mouse data','Multiselect','on');

%%
clear HPC RSC
if ischar(fnm)
    fnm={fnm};
end
lg=[];
h=1; r=1;
for i=1:length(fnm)
    load([pth fnm{1,i}])
    sp=split(fnm{1,i},'.');
    switch sp{1,1}(1,end)
        case 'H'
    HPC{h}=NF_map.list;
    h=h+1;
        case 'R'
          RSC{r}=NF_map.list;
          r=r+1;
    end
end

%% 'Fraction of TXN Cell (%)'
clear SEM
    figure(1)
   p=[];

for i=1:size(HPC,2)
    fracTXN{1}(i,1)=sum(HPC{i}(:,3)==2)/sum(HPC{i}(:,3)~=3);
    numTXN{1}(i,1)=sum(HPC{i}(:,4)==2)/sum(HPC{i}(:,3)==2);   
end
SEM(1,1)=std(fracTXN{1},0,1,'omitnan')/sqrt(size(HPC,2));
SEM(2,1)=std(numTXN{1},0,1,'omitnan')/sqrt(size(HPC,2));

for i=1:size(RSC,2)
    fracTXN{2}(i,1)=sum(RSC{i}(:,3)==2)/sum(RSC{i}(:,3)~=3);
    numTXN{2}(i,1)=sum(RSC{i}(:,4)==2)/sum(RSC{i}(:,3)==2);   
end
SEM(1,2)=std(fracTXN{2},0,1,'omitnan')/sqrt(size(RSC,2));
SEM(2,2)=std(numTXN{2},0,1,'omitnan')/sqrt(size(RSC,2));

bar([1:1:2],[mean(fracTXN{1},1) mean(fracTXN{2},1)]*100,'FaceColor',[0.7 0.7 0.7],'BarWidth',0.8)%'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:1:2],'XTickLabel',...
    {'HPC','RSC'});
   ylabel('Fraction of TXN Cell (%)')
   xlim([0.5 2.5])
   hold all
   errorbar([1:1:2],[mean(fracTXN{1},1) mean(fracTXN{2},1)]*100,SEM(1,:)*100,'k','linewidth',4,'linestyle','none')
xtickangle(45);


figure(2)
bar([1:1:2],[mean(numTXN{1},1) mean(numTXN{2},1)]*100,'FaceColor',[0.7 0.7 0.7],'BarWidth',0.8)%'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:1:2],'XTickLabel',...
    {'HPC','RSC'});
   ylabel('Fraction of cells with 2 TXN site (%)')
   xlim([0.5 2.5])
   hold all
   errorbar([1:1:2],[mean(numTXN{1},1) mean(numTXN{2},1)]*100,SEM(2,:)*100,'k','linewidth',4,'linestyle','none')
xtickangle(45);

