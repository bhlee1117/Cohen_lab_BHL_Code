clear
[fnma ptha]=uigetfile('*.mat','Multiselect','On');
load('C:\Users\Administrator\Dropbox\In vivo imaging\Analysis\Blind_Answer_total.mat');
mother_pth='C:\Users\Administrator\Dropbox\In vivo imaging\Analysis\';
%%
clear data
for i=1:size(fnma,2) %Mouse number
    sp=split(fnma{1,i},'.');
    list{1}=dir([mother_pth '\' sp{1,1} '\' 'Blinder1\']);
    list{2}=dir([mother_pth '\' sp{1,1} '\' 'Blinder2\']);
    load([ptha fnma{1,i}])
    
    for j=1:size(BlindAnswer,1)
        if ~isempty(strfind(sp{1,1},char(BlindAnswer(j,1)))) && ~isempty(strfind(char(BlindAnswer(j,1)),sp{1,1}))
            ref=j; %row number that matches to the mouse number
        end
    end
            
    for j=3:size(list{1},1) %image number
    data{i,j-2}(:,3)=Final.Data(:,j+1);
    
    load([mother_pth '\' sp{1,1} '\' 'Blinder1\' list{1}(j).name])
    data{i,find(double(BlindAnswer(ref,2:end))==j-2)}(:,1)=NF_map.list(:,4);
    load([mother_pth '\' sp{1,1} '\' 'Blinder2\' list{1}(j).name])
    data{i,find(double(BlindAnswer(ref,2:end))==j-2)}(:,2)=NF_map.list(:,4);
 
    end
end
%%
clear accuracy
for i=1:size(data,1)
    for j=2:size(data,2)
        if ~isempty(data{i,j}) && size(data{i,j},2)>2 && isempty(find(data{i,j}==0))
    accuracy{1}(i,j)=[sum(data{i,j}(:,1)==data{i,j}(:,3))/size(data{i,j},1)];
    accuracy{2}(i,j)=sum(data{i,j}(:,2)==data{i,j}(:,3))/size(data{i,j},1);
    accuracy{3}(i,j)=[sum(data{i,j}(:,1)==data{i,j}(:,2))/size(data{i,j},1)];                    
 accuracy{4}(i,j)=sum((data{i,j}(:,1)==2)&(data{i,j}(:,3)==2))/sum(data{i,j}(:,3)==2);
 accuracy{5}(i,j)=sum((data{i,j}(:,2)==2)&(data{i,j}(:,3)==2))/sum(data{i,j}(:,3)==2);
        else


    end
    end
end
%%
clear mean_acc
g=1;
for i=[1:7 9:23]
    for j=1:5
    accuracy{j}(accuracy{j}==0)=NaN;
    mean_acc(g,j)=mean(accuracy{j}(i,:),'omitnan');
    end
    g=g+1;
end

for i=1:size(mean_acc,1)
    plot([1:size(mean_acc,2)]+randn(1,5)*0.05,mean_acc(i,:),'color',[0.3 0.3 0.3],'linestyle','none','marker','.')
    hold all
end
plot([1:size(mean_acc,2)],mean(mean_acc,1),'ro')
errorbar([1:size(mean_acc,2)],mean(mean_acc,1),std(mean_acc),'color',[1 0 0],'linestyle','none','marker','.')
set(gca,'FontSize',8,'LineWidth',1,'XTick',[1:size(mean_acc,2)],'XTickLabel',...
    {'Blinder1 & Final','Blinder2 & Final','Blinder1 & Blinder2','Blinder1 & Final (only TXN)','Blinder2 & Final (only TXN)'});
xlim([0.5 size(mean_acc,2)+0.5])
ylim([0.2 1])
xtickangle(45)