clear
clc;
cd '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/FromBackup/PP72_PlaceCellResults';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:P23');

% [~, ~, NeuronsToUse]=xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
%     'PlaceCellData_Arrangement.xlsx'], 'Sheet1', 'L8:M46');
%
% NeuronsToUse=cellfun(@(x) (str2num(num2str(x))),NeuronsToUse,'UniformOutput',false);
ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
fpath=raw(:,1)';
StructureData=raw(:,10);
EndFrame=cell2mat(raw(:,13));
ifmotionReject=cell2mat(raw(:,14));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
%%
f=18;
load([fpath{f} '/Result.mat'])
Sp_time=find(max(Result.SpClass(1:2,:),[],1));
[~, Sp_type]=max(Result.SpClass(1:2,Sp_time),[],1);
Spike_period=max([Result.SpClass(1,:); Result.CStrace]);
for s=1:2
precedeing_spike{s}=[];    
    for sp=find(Sp_type(2:end)==s)
Spike_period2=bwlabel(Spike_period);
Spike_period2(Spike_period2==Spike_period2(Sp_time(sp+1)))=0;
Spike_period2=Spike_period2~=0;
last_spike=max(find(Spike_period2(1:Sp_time(sp+1)-1)));
precedeing_spike{s}=[precedeing_spike{s}; [Sp_time(sp+1)-last_spike Sp_type(sp)]];
    end
end

figure(3); clf; clear h
tiledlayout(2,2); cmap=distinguishable_colors(2);
bin_edge=10.^[0:0.2:5];
title_str={'Preceding spike before SS','Preceding spike before CS'};
for s=1:2
nexttile([1 1]);
for ss=1:2
ind=find(precedeing_spike{s}(:,2)==ss);
h=histogram(precedeing_spike{s}(ind,1),bin_edge,'Normalization','probability'); hold all
end
set(gca,'XScale','log')
xlabel('Preceding silence (ms)')
ylabel('Probability')
title(title_str{s})
legend({'SS preceding','CS preceding'})
end

nexttile([1 1])
scatter(1+randn(size(precedeing_spike{1},1),1)*0.1,precedeing_spike{1}(:,1),10,cmap(1,:)); hold all
scatter(2+randn(size(precedeing_spike{2},1),1)*0.1,precedeing_spike{2}(:,1),10,cmap(2,:))
errorbar([1 2],[mean(precedeing_spike{1}(:,1)) mean(precedeing_spike{2}(:,1))],[std(precedeing_spike{1}(:,1)) std(precedeing_spike{2}(:,1))],'color','k')
set(gca,'yScale','log','XTick',[1 2],'XTickLabel',{'SS','CS'})
ylabel('Preceding silence (ms)')

nexttile([1 1])
ratio_CS=[];
for b=1:length(bin_edge)-1
    for s=1:2
    data_inBin=cellfun(@(x) find(x(:,1)>=bin_edge(b) & x(:,1)<bin_edge(b+1) & x(:,2)==s),precedeing_spike,'UniformOutput',false);
    ratio_CS(b,s)=length(data_inBin{2})/(length(data_inBin{1})+length(data_inBin{2}));
    end
end
plot(mean([bin_edge(1:end-1); bin_edge(2:end)]),ratio_CS(:,1)); hold all
plot(mean([bin_edge(1:end-1); bin_edge(2:end)]),ratio_CS(:,2));
set(gca,'XScale','log')
xlabel('Preceding silence (ms)')
ylabel('Probability to show CS next')
legend({'SS','CS'})
%%
load('/Volumes/BHL18TB_D2/YQ601_PlaceCellResults/PF_Result_20240111.mat')
load('/Volumes/BHL18TB_D2/Arranged_Data/NaV_InactivationResult/PlaceCell_List.mat')
%%
for i=1:size(PlaceCell_List,1)
spike=PC_Result{PlaceCell_List(i,1)}.spike(PlaceCell_List(i,2),:);
[~, CS_time]=unique(spike.*bwlabel(PC_Result{PlaceCell_List(i,1)}.CS_trace(PlaceCell_List(i,2),:)));
CS_time=CS_time(2:end);
SpClass=zeros(2,size(PC_Result{PlaceCell_List(i,1)}.traces,2));
SpClass=spike.*(1-PC_Result{PlaceCell_List(i,1)}.CS_trace(PlaceCell_List(i,2),:));
SpClass(2,CS_time)=1;
CStrace=PC_Result{PlaceCell_List(i,1)}.CS_trace(PlaceCell_List(i,2),:);

Sp_time=find(max(SpClass(1:2,:),[],1));
[~, Sp_type]=max(SpClass(1:2,Sp_time),[],1);
Spike_period=max([SpClass(1,:); CStrace]);

for s=1:2
precedeing_spike{i,s}=[];    
    for sp=find(Sp_type(2:end)==s)
Spike_period2=bwlabel(Spike_period);
Spike_period2(Spike_period2==Spike_period2(Sp_time(sp+1)))=0;
Spike_period2=Spike_period2~=0;
last_spike=max(find(Spike_period2(1:Sp_time(sp+1)-1)));
precedeing_spike{i,s}=[precedeing_spike{i,s}; [Sp_time(sp+1)-last_spike Sp_type(sp)]];
    end
end
end
%%
bin_edge=10.^[0:0.4:5];
histsilence=[]; histsilence_cat=[];
for i=1:size(PlaceCell_List,1)
    for s=1:2
        for ss=1:2
            ind=find(precedeing_spike{i,s}(:,2)==ss);
            histsilence((s-1)*2+ss,:,i)=histcounts(precedeing_spike{i,s}(ind,1),bin_edge,'Normalization','probability');            
        end
        histsilence_cat(s,:,i)=histcounts(precedeing_spike{i,s}(:,1),bin_edge,'Normalization','probability');
    end
end
%%

Meanhistsilence=mean(histsilence,3);
Stdhistsilence=std(histsilence,0,3,'omitnan');

figure(3); clf; 
tiledlayout(2,2); cmap=distinguishable_colors(2);
xtick_bin=mean([bin_edge(1:end-1); bin_edge(2:end)]);
title_str={'Preceding spike before SS','Preceding spike before CS'};
for s=1:2
nexttile([1 1]); clear l
for ss=1:2
l(ss)=errorbar_shade(xtick_bin,Meanhistsilence((s-1)*2+ss,:),Stdhistsilence((s-1)*2+ss,:),cmap(ss,:)); hold all
end
set(gca,'XScale','log')
xlabel('Preceding silence (ms)')
ylabel('Probability')
title(title_str{s})
legend(l,{'SS preceding','CS preceding'})
end

Meanhistsilence_cat=mean(histsilence_cat,3);
Stdhistsilence_cat=std(histsilence_cat,0,3,'omitnan');
nexttile([1 1]); clear l
for s=1:2
l(s)=errorbar_shade(xtick_bin,Meanhistsilence_cat(s,:),Stdhistsilence_cat(s,:),cmap(s,:)); hold all
end
legend(l,{'SS','CS'})
set(gca,'XScale','log')
xlabel('Preceding silence (ms)')
ylabel('Probability')

nexttile([1 1])
ratio_CS=[];
for b=1:length(bin_edge)-1
    for s=1:2
    data_inBin=cellfun(@(x) find(x(:,1)>=bin_edge(b) & x(:,1)<bin_edge(b+1) & x(:,2)==s),precedeing_spike,'UniformOutput',false);
    for i=1:size(PlaceCell_List,1)
    ratio_CS(b,s,i)=length(data_inBin{i,2})/(length(data_inBin{i,1})+length(data_inBin{i,2}));
    end
    end
end
Meanratio_CS=mean(ratio_CS,3,'omitnan'); Stdratio_CS=std(ratio_CS,0,3,'omitnan');
Meanratio_CS(isnan(Meanratio_CS))=0; Stdratio_CS(isnan(Stdratio_CS))=0;
%plot(xtick_bin,squeeze(ratio_CS(:,1,:))','color',[0.3 0.7 1]); hold all
%plot(xtick_bin,squeeze(ratio_CS(:,2,:))','color',[1 0.67 0.85]); hold all
l(1)=errorbar_shade(xtick_bin,Meanratio_CS(:,1),Stdratio_CS(:,1),cmap(1,:));
l(2)=errorbar_shade(xtick_bin,Meanratio_CS(:,2),Stdratio_CS(:,2),cmap(2,:));
set(gca,'XScale','log')
xlabel('Preceding silence (ms)')
ylabel('Probability to show CS next')
legend(l,{'SS','CS'})


