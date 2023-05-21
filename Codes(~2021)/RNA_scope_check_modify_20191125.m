%% Load mouse data

clear
[fnm pth]=uigetfile('.mat','Select the analyzed mouse data','Multiselect','on');
lg{1}=[{'CFC-2-1'},{'CFC-2-2'},{'CFC-2-3'},{'CFC-3-1'},{'CFC-3-2'},{'CFC-3-3'}];
lg{2}=[{'HC-1-1'},{'HC-1-2'},{'HC-1-3'},{'HC-2-1'},{'HC-2-2'},{'HC-2-3'},{'HC-3-2'},{'HC-3-3'}];
%%
clear CFC HC
if ischar(fnm)
    fnm={fnm};
end
lg=[];
h=1; r=1;
load(['H:\Image_data\KIST_RNAscope\20191101_SNU-RNAscope_Arc-Fos-Egr1\Analysis\Classification\CFC_modifiy.mat'])
load(['H:\Image_data\KIST_RNAscope\20191101_SNU-RNAscope_Arc-Fos-Egr1\Analysis\Classification\HC_modifiy.mat'])
for i=1:length(fnm)
    load([pth fnm{1,i}])
    
    sp=split(fnm{1,i},'_');
    if i<7
        CFC_mod{1,i}(CFC_mod{1,i}==1)=2;
        CFC_mod{1,i}(CFC_mod{1,i}==0)=1;
for j=1:size(CFC_mod{1,i},1)    
RNAscope.result(CFC_mod{1,i}(j,1),3:5)=CFC_mod{1,i}(j,2:4);
end
    end
    
    if i>=7
        HC_mod{1,i-6}(HC_mod{1,i-6}==1)=2;
        HC_mod{1,i-6}(HC_mod{1,i-6}==0)=1;
for j=1:size(HC_mod{1,i-6},1)    
RNAscope.result(HC_mod{1,i-6}(j,1),3:5)=HC_mod{1,i-6}(j,2:4);
end
    end
    
    switch sp{3,1}(1)
        case 'C'
    
    CFC{1,h}=RNAscope.result(find(sum(RNAscope.result(:,3:5)==3,2)==0),:); % list without NA cells
    CFC{2,h}=RNAscope.result(find(sum(RNAscope.result(:,3:5)==3,2)~=0),:); % list of NA cells
    CFC{3,h}=RNAscope.max; % list with NA cells
    CFC{4,h}=find(sum(RNAscope.result(:,3:5)==3,2)==0);
    CFC{5,h}=find(sum(RNAscope.result(:,3:5)==3,2)~=0);
    h=h+1;
        case 'H'
          HC{1,r}=RNAscope.result(find(sum(RNAscope.result(:,3:5)==3,2)==0),:);
          HC{2,r}=RNAscope.result; % list with NA cells
          HC{3,r}=RNAscope.max;
          HC{4,r}=find(sum(RNAscope.result(:,3:5)==3,2)==0);
          HC{5,r}=find(sum(RNAscope.result(:,3:5)==3,2)~=0);
          r=r+1;
    end
end

%%
for i=3%:size(CFC,2)
       handles.fig = figure(1);
set(gcf, 'Position',  [50, 30, 5000, 2000])
color=[0.7 .7 .7;0 1 0;2 0 0;0 0.9 0.9];
    RGB_im=zeros(size(CFC{3,i},1),size(CFC{3,i},2),3);
    for col=1:size(color,1)
        for rgb=1:size(color,2)
            CFC{3,i}=double(CFC{3,i});
    RGB_im(:,:,rgb)=RGB_im(:,:,rgb)+CFC{3,i}(:,:,col)/max(max(CFC{3,i}(:,:,col)))*color(col,rgb);
        end
    end
    imagesc(RGB_im)
     axis equal tight off
%     imagesc(CFC{3,i}(:,:,1))
%     axis equal tight
%     hold all
%     for j=1:size(CFC{1,i},1)
%     plot(CFC{1,i}(j,1),CFC{1,i}(j,2),'wo','markersize',20)
%     text(CFC{1,i}(j,1)-10,CFC{1,i}(j,2),num2str(CFC{4,i}(j,1)),'color','w')
%     end
%     for j=1:size(CFC{2,i},1)
%     plot(CFC{2,i}(:,1),CFC{2,i}(:,2),'ro','markersize',20)
%     text(CFC{2,i}(j,1)-10,CFC{2,i}(j,2),num2str(CFC{5,i}(j,1)),'color','w')
%     end
    
    F = getframe(gcf);
[X, Map] = frame2im(F);
close(figure(1))

    
        handles.fig = figure(2);
set(gcf, 'Position',[50, 30, 5000, 2000])
    imagesc(RGB_im)
    axis equal tight off
    hold all
    for j=1:size(CFC{1,i},1)
    plot(CFC{1,i}(j,1),CFC{1,i}(j,2),'wo','markersize',20)
    text(CFC{1,i}(j,1)-10,CFC{1,i}(j,2)-40,num2str(CFC{4,i}(j,1)),'color','w')
    end
    color_mark=[0 1 0 ;1 0 0 ;0 1 1];
    for c=3:5
    for j=find(CFC{1,i}(:,c)==2)
    text(CFC{1,i}(j,1)+7*(c-4),CFC{1,i}(j,2),'o','color',color_mark(c-2,:))
    
    end
    end
    FF = getframe(gcf);
[XX, Map2] = frame2im(FF);
    close(figure(2))
    my_image_compare_tool(X,XX);
    
    
end
%%
for i=1:size(HC,2)
end

%% 'Fraction of TXN Cell (%)'
clear SEM
    figure(1)
   p=[];

for i=1:size(CFC,2)
    for c=1:3
    fracTXN{1}(i,c)=sum(CFC{1,i}(:,2+c)==2)/sum(CFC{1,i}(:,2+c)~=3); %{1} for CFC
    end
end
SEM(1,:)=std(fracTXN{1},0,1,'omitnan')/sqrt(size(CFC,2));

for i=1:size(HC,2)
    for c=1:3
        fracTXN{2}(i,c)=sum(HC{1,i}(:,2+c)==2)/sum(HC{1,i}(:,2+c)~=3); %{2} for HC
    end
end
SEM(2,:)=std(fracTXN{2},0,1,'omitnan')/sqrt(size(HC,2));


bar([1 2],[mean(fracTXN{1},1) ; mean(fracTXN{2},1)]*100,'FaceColor',[0.7 0.7 0.7],'BarWidth',0.8)%'FaceColor',[0.6 0.6 0.6],'BarWidth',0.8)
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1 2],'XTickLabel',...
    {'CFC','HC'});
   ylabel('Fraction of TXN Cell (%)')
   xlim([0.5 2.5])
   hold all
   errorbar([0.77 1 1.23 1.77 2 2.23],[mean(fracTXN{1},1) mean(fracTXN{2},1)]*100,[SEM(1,:)*100 SEM(2,:)*100],'k','linewidth',4,'linestyle','none')
xtickangle(45);
%%
color=[0 1 0;1 0 0;0 1 1;0 0.45 0.74;0.85 0.33 0.1;0.93 0.69 0.13;0.49 0.18 0.56;0.47 0.67 0.19;0.64 0 0.8];
gr=0.23;
for i=1:2 %% plot by line
    for j=1:size(fracTXN{i},1)
        plot([i-gr i i+gr],fracTXN{i}(j,1:end)*100,'.','markersize',15,'linestyle','-','color',color(j,:))
        hold all
        text(i+gr+0.1,fracTXN{i}(j,end)*100,lg{i}{j},'color',color(j,:))
    end
end
errorbar([0.77 1 1.23],[mean(fracTXN{1},1)]*100,[SEM(1,:)*100],'k','linewidth',2)
errorbar([1.77 2 2.23],[mean(fracTXN{2},1)]*100,[SEM(2,:)*100],'k','linewidth',2)
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1-gr 1 1+gr 2-gr 2 2+gr],'XTickLabel',...
    {'Arc CFC','C-fos CFC','Egr1 CFC','Arc HC','C-fos HC','Egr1 HC'});
   ylabel('Fraction of TXN Cell (%)')
   xlim([0.5 3])
   xtickangle(45);
%%

clear overlap over_name overlapTXN SEM data
Window_size=200;
exhib=[3 4 5];
for i=[1] %CFC=1 HC=2
%HC CFC Ret1
  handles.fig = figure(1);
set(gcf, 'Position',  [50, 30, Window_size*(size(CFC,2)+2), 300])

for j=1:size(CFC,2)
  
handles.axes1 = axes('Units','pixels','Position',[Window_size*(j-1)+50 50 Window_size Window_size]);   
title(lg{i}{j})
tmp=venn_data(CFC{1,j}(:,exhib))/sum(sum(CFC{1,j}==2,2)>0);%size(CFC{1,j},1);
data{i}(j,:)=tmp;
if j~=6
    [a b]=venn(data{i}(j,:));
end
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
clear cond_prob
exhib=[3:5];
overlap=[2 1; 3 1;3 2;4 1;4 2;4 3];
%Arc:1, C-fos:2, Egr-1:3, All:4
Dat{1}=CFC;
Dat{2}=HC;
overlap=overlap+2;
for i=1:2
    for k=1:size(Dat{i},2) %mouse
    for j=1:size(overlap,1) %conditions
   
        if overlap(j,1)==6
        cond_prob{i}(j,k)=sum(sum(Dat{i}{1,k}(:,3:5)==2,2)==3)/sum(Dat{i}{1,k}(:,overlap(j,2))==2);  
        else
        cond_prob{i}(j,k)=sum(Dat{i}{1,k}(:,overlap(j,1))==2&Dat{i}{1,k}(:,overlap(j,2))==2)/sum(Dat{i}{1,k}(:,overlap(j,2))==2);  
        end
    end
    end
     for k=1:size(Dat{i},2) %mouse
        if i==1
            plot([1:size(overlap,1)],cond_prob{i}(:,k),'.','markersize',15,'color',[0.8 0.7 0.7],'linestyle','-')
        else
            plot([1:size(overlap,1)],cond_prob{i}(:,k),'.','markersize',15,'color',[0.7 0.8 0.7],'linestyle','-')
        end
        hold all
   set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:size(overlap,1)],'XTickLabel',...
    {'C-fos & Arc | Arc', 'Arc & Egr-1 | Arc', 'Egr-1 & C-fos | C-fos', 'All | Arc', 'All | C-fos', 'All | Egr-1'});
   ylabel('Conditional probability (%)')
   xtickangle(45);
   xlim([0.5 6.5])
     end 
end  
    for i=1:2
    if i==1 
         plot([1:size(overlap,1)],mean(cond_prob{i},2),'.','markersize',15,'color',[1 0 0],'linestyle','-')
     else
          plot([1:size(overlap,1)],mean(cond_prob{i},2),'.','markersize',15,'color',[0 1 0],'linestyle','-')
     end
end