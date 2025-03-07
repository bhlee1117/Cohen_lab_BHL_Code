clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:Q31');

ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
basal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
apical_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false);
fpath=raw(:,1)';
StructureData=raw(:,10);
BadROI=cellfun(@(x) (str2num(num2str(x))),raw(:,15),'UniformOutput',false);
EndFrame=cell2mat(raw(:,13));
ifmotionReject=cell2mat(raw(:,14));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
set(0,'DefaultFigureWindowStyle','docked')

%%
f=13; load(fullfile(fpath{f},'PC_Result.mat'),'Result')
cd(fpath{f})
rois={basal_ROI{f},apical_ROI{f},ref_ROI{f}};
nROI=size(Result.normTraces,1);
nTau_bAP=[-20:20];
nTau={[-120:50],[-120:50],[-30:20],[-120:50]}; %SS, CS, dSP, Brst
StimOn_Lap=unique(Result.VR(8,:).*double(Result.Blue>0)); StimOn_Lap=StimOn_Lap(2:end);
StimOff_Lap=unique(Result.VR(8,:).*double(Result.VR(2,:)==2)); StimOff_Lap=setdiff(StimOff_Lap,[0 StimOn_Lap]);
filterfreq=[];
% Isolated Somatic spike
som_spike=find(Result.spike(1,:)); 
ss_time=find(Result.SpClass(1,:)); brst=bwlabel((ss_time(2:end)-ss_time(1:end-1))<15);
SpClass=Result.SpClass; BS_trace=zeros(1,size(Result.traces,2));
for b=1:max(brst)
    bwn=find(brst==b);
    SpClass(1,ss_time([bwn bwn(end)+1]))=0;
    SpClass(4,ss_time([bwn(1)]))=1;
    BS_trace(1,[ss_time(bwn): ss_time(bwn(end)+1)])=b;
end
tr_ref=Result.normTraces(ref_ROI{f},:);
tr_sub=mean(tr_ref,1)-movprc(mean(tr_ref,1),200,20,2);
tr_sub=get_subthreshold(tr_sub,Result.spike(1,:),5,10);
[trans tr_trace]=detect_transient2(tr_sub,[5 1.5],Result.spike(1,:),15);

dist_order=Result.dist_order;

bAP_ref=[];
for s=som_spike
    isnearby=sum(ismember(s+nTau_bAP,som_spike))>1;
    isnearbyCS=sum(ismember(s+nTau_bAP,find(Result.CStrace)))>1;
    ispartCS=tr_trace(s)>0;
    if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
        bAP_ref=[bAP_ref s];
    end
end
sp_na=sum((bAP_ref'+nTau_bAP)<0 | (bAP_ref'+nTau_bAP)>size(Result.traces,2),2)==0;
bAP_ref=bAP_ref(sp_na);

% Isolated Somatic spike
bAP_s=[];
for s=som_spike
    isnearby=sum(ismember(s+nTau{1},som_spike))>1;
    isnearbyCS=sum(ismember(s+nTau{1},find(Result.CStrace)))>1;
    ispartCS=tr_trace(s)>0;
    if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
        bAP_s=[bAP_s s];
    end
end
sp_na=sum((bAP_s'+nTau{1})<0 | (bAP_s'+nTau{1})>size(Result.traces,2),2)==0;
bAP_s=bAP_s(sp_na);

dSP_s=[];
dSP_list=find(Result.SpClass(3,:));
isApical=sum(reshape(Result.spike(rois{2},dSP_list'+[-1:1]),length(rois{2}),[],3),[1 3])>0;
isBasal=sum(reshape(Result.spike(rois{1},dSP_list'+[-1:1]),length(rois{1}),[],3),[1 3])>0;
isBlueoff=Result.Blue(dSP_list)==0;
isSoma=sum(reshape(Result.spike(rois{3},dSP_list'+[-1:1]),length(rois{3}),[],3),[1 3])>0;
dSP_list=dSP_list(isBlueoff & ~isBasal);
for s=dSP_list
    isnearby=sum(ismember(s+nTau{3},som_spike))>1;
    isnearbyCS=sum(ismember(s+nTau{3},find(Result.CStrace)))>1;
    ispartCS=tr_trace(s)>0;
    if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
        dSP_s=[dSP_s s];
    end
end
sp_na=sum((dSP_s'+nTau{3})<0 | (dSP_s'+nTau{3})>size(Result.traces,2),2)==0;
dSP_s=dSP_s(sp_na);

% Isolated Complex spike
CS_s=[];
C_list=find(Result.SpClass(2,:));
CS_label=max(Result.spike,[],1).*bwlabel(Result.CStrace);
CS_list=[];
for g=1:length(C_list)
    s=C_list(g);
    s_tmp=Result.spike(1,:);
    s_tmp(find(CS_label==g))=0;
    isnearbyCS=max(bwlabel(Result.CStrace(s+nTau{2})))>1;
    isnearbyS=sum(ismember(s+nTau{2},find(s_tmp)))>0;
    if ~isnearbyCS & ~isnearbyS
        CS_s=[CS_s s];
        CS_list=[CS_list g];
    end
end
sp_na=sum((CS_s'+nTau{2})<0 | (CS_s'+nTau{2})>size(Result.traces,2),2)==0;
CS_s=CS_s(sp_na);

BS_s=[];
B_list=find(SpClass(4,:));
for g=1:length(B_list)
    s=B_list(g);
    
    isnearby=sum(ismember(s+nTau{4},find(BS_trace~=BS_trace(s) & Result.spike(1,:)==1)))>0;
    isnearbyCS=sum(ismember(s+nTau{4},find(Result.CStrace)))>0;
    if ~isnearbyCS & ~isnearby
        BS_s=[BS_s s];
    end
end
sp_na=sum((BS_s'+nTau{4})<0 | (BS_s'+nTau{4})>size(Result.traces,2),2)==0;
BS_s=BS_s(sp_na);


STA_SS=squeeze(mean(reshape(Result.normTraces(:,bAP_s'+nTau_bAP),nROI,[],length(nTau_bAP)),2));
STA_SS= STA_SS - mean(STA_SS(:,1:10),2);
%STA_SS= STA_SS - prctile(STA_SS,30,2);
F_ref=mean(STA_SS(:,-nTau_bAP(1)+[7:11]),2);
%F_ref=(tovec(imgaussfilt(Result.ref_im,1))'*tovec(Result.ftprnt)/Result.SpikeHeight_fit(1))';
%F_ref=median(-Result.traces_bvMask(:,1:2000),2);

SilentPeriod=ones(1,size(Result.traces,2));
sp_time=find(max(Result.spike,[],1))';
sp_na=sum((find(max(Result.spike,[],1))'+[-10:150])<0 | (find(max(Result.spike,[],1))'+[-10:150])>size(Result.traces,2),2)==0;
SilentPeriod(sp_time(sp_na)+[-10:150])=NaN;
t_fit=find(~isnan(SilentPeriod) & Result.Blue==0);

NormalizedTrace=(Result.normTraces)./F_ref;
if ~isempty(filterfreq)
    freq_lowhigh=filterfreq/(1000/2);
    [b, a] = butter(4, freq_lowhigh, 'stop');
    for n=1:size(NormalizedTrace,1)
        NormalizedTrace(n,:) = filtfilt(b, a, NormalizedTrace(n,:)');
    end
end
% lwpass=NaN(size(NormalizedTrace));
% lwpass(:,t_fit)=NormalizedTrace(:,t_fit);
% lwpass=movmedian(lwpass,30000,2,'omitnan');
% NormalizedTrace=NormalizedTrace-lwpass;
NormalizedTrace_dirt=NormalizedTrace;
%NormalizedTrace_dirt(:,17000:18000)=NaN;
NormalizedTrace_dirt(:,Result.motionReject)=NaN;
NormalizedTrace_ch=cellfun(@(x) x./F_ref,Result.norm_trace_check,'UniformOutput',false);
NormalizedTrace_ch{1}(:,Result.motionReject)=NaN; NormalizedTrace_ch{2}(:,Result.motionReject)=NaN;

if isfield(Result,'dirtTrace')
NormalizedTrace_dirt(Result.dirtTrace>0)=NaN;
NormalizedTrace_ch{1}(Result.dirtTrace>0)=NaN; NormalizedTrace_ch{2}(Result.dirtTrace>0)=NaN;
end

STA_CSmat=reshape(NormalizedTrace_dirt(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));
STA_SSmat=reshape(NormalizedTrace_dirt(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1}));
STA_dSPmat=reshape(NormalizedTrace_dirt(:,dSP_s'+nTau{3}),nROI,[],length(nTau{3}));
STA_BSmat=reshape(NormalizedTrace_dirt(:,BS_s'+nTau{4}),nROI,[],length(nTau{4}));

STA_CSmat_ch=cellfun(@(x) reshape(x(:,CS_s'+nTau{2}),nROI,[],length(nTau{2})),NormalizedTrace_ch,'UniformOutput',false);
STA_SSmat_ch=cellfun(@(x) reshape(x(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1})),NormalizedTrace_ch,'UniformOutput',false);
STA_dSPmat_ch=cellfun(@(x) reshape(x(:,dSP_s'+nTau{3}),nROI,[],length(nTau{3})),NormalizedTrace_ch,'UniformOutput',false);
STA_BSmat_ch=cellfun(@(x) reshape(x(:,BS_s'+nTau{4}),nROI,[],length(nTau{4})),NormalizedTrace_ch,'UniformOutput',false);
%%
figure(4); clf; cmap=distinguishable_colors(6); ax1=[];
noi=setdiff([1:nROI],[BadROI{f}]);
noi_dist=ismember(Result.dist_order,noi);
%rois={basal_ROI{f},apical_ROI{f}};
cax=[-1 10];%*0.01;
tiledlayout(2,4);
ax3=nexttile([1 1]);
imagesc(nTau{1},[1:sum(noi_dist)],squeeze(mean(STA_SSmat(dist_order(noi_dist),:,:)-median(STA_SSmat(dist_order(noi_dist),:,1:-nTau{1}(1)),3,'omitnan'),2,'omitnan')),cax)
title('Simple spike')
ax4=nexttile([1 1]);
imagesc(nTau{2},[1:sum(noi_dist)],squeeze(mean(STA_CSmat(dist_order(noi_dist),:,:)-median(STA_CSmat(dist_order(noi_dist),:,1:-nTau{2}(1)),3,'omitnan'),2,'omitnan')),cax); %shading interp;
title('Complex spike')
ax5=nexttile([1 1]);
imagesc(nTau{3},[1:sum(noi_dist)],squeeze(mean(STA_dSPmat(dist_order(noi_dist),:,:)-median(STA_dSPmat(dist_order(noi_dist),:,1:-nTau{3}(1)),3,'omitnan'),2,'omitnan')),cax)
title('Dendritic spike')
ax6=nexttile([1 1]);
imagesc(nTau{4},[1:sum(noi_dist)],squeeze(mean(STA_BSmat(dist_order(noi_dist),:,:)-median(STA_BSmat(dist_order(noi_dist),:,1:-nTau{4}(1)),3,'omitnan'),2,'omitnan')),cax)
title('Burst spike')
colormap(turbo); linkaxes([ax3 ax4 ax5 ax6],'xy')

ax1=[ax1 nexttile([1 1])];
l=plot(nTau{1},squeeze(mean(STA_SSmat(dist_order(noi_dist),:,:)-median(STA_SSmat(dist_order(noi_dist),:,1:-nTau{1}(1)),3,'omitnan'),2,'omitnan'))); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
ax1=[ax1 nexttile([1 1])];
l=plot(nTau{2},squeeze(mean(STA_CSmat(dist_order(noi_dist),:,:)-median(STA_CSmat(dist_order(noi_dist),:,1:-nTau{2}(1)),3,'omitnan'),2,'omitnan'))); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
ax1=[ax1 nexttile([1 1])];
l=plot(nTau{3},squeeze(mean(STA_dSPmat(dist_order(noi_dist),:,:)-median(STA_dSPmat(dist_order(noi_dist),:,1:-nTau{3}(1)),3,'omitnan'),2,'omitnan'))); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
ax1=[ax1 nexttile([1 1])];
l=plot(nTau{4},squeeze(mean(STA_BSmat(dist_order(noi_dist),:,:)-median(STA_BSmat(dist_order(noi_dist),:,1:-nTau{4}(1)),3,'omitnan'),2,'omitnan'))); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
linkaxes([ax1],'xy')

figure(6); clf;
blueOn_sp=mean(SpClass(:,Result.Blue>0 & Result.motionReject==0),2,'omitnan')*1000;
blueOff_sp=mean(SpClass(:,Result.Blue==0 & Result.motionReject==0),2,'omitnan')*1000;
b=bar([1:4],[blueOn_sp blueOff_sp]','BarLayout','grouped','FaceColor','flat');
b(1).CData = repmat([0 0.6 1],4,1); b(2).CData = repmat([0.3 0.3 0.3],4,1);
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'SS','CS','dSP','Brst'},'yscale','log')
ylabel('Firing rate (Hz)');
legend({'Blue On','Blue Off'})

figure(7); clf;
tiledlayout(1,4);
EdgeBin=[-4:0.4:4];
SS_s_off=find(Result.Blue(bAP_s)==0);
CS_s_off=find(Result.Blue(CS_s)==0);
BS_s_off=find(Result.Blue(BS_s)==0);
Prespike_ampBS=mean(STA_BSmat(Result.dist_order(noi_dist),BS_s_off,-nTau{4}(1)+[-10:-1]),3,'omitnan');
Prespike_ampCS=mean(STA_CSmat(Result.dist_order(noi_dist),CS_s_off,-nTau{2}(1)+[-10:-1]),3,'omitnan');
Prespike_ampSS=mean(STA_SSmat(Result.dist_order(noi_dist),SS_s_off,-nTau{1}(1)+[-10:-1]),3,'omitnan');

somROI=find(Result.dist_order(noi_dist)==1);
nexttile(1,[1 1])
showdat={mean(Prespike_ampSS(somROI+[-1:1],:),1,'omitnan')',mean(Prespike_ampCS(somROI+[-1:1],:),1,'omitnan')',mean(Prespike_ampBS(somROI+[-1:1],:),1,'omitnan')'};
plot_errorbar3([1:3],showdat,'ranksum',[-4 5],'Pre-spike AUC, peri-soma',{'SS','CS','BS'},distinguishable_colors(3));

nexttile(2,[1 1])
showdat={mean(Prespike_ampSS([1:3],:),1,'omitnan')',mean(Prespike_ampCS([1:3],:),1,'omitnan')',mean(Prespike_ampBS([1:3],:),1,'omitnan')'};
plot_errorbar3([1:3],showdat,'ranksum',[-4 5],'Pre-spike AUC, basal dendrite',{'SS','CS','BS'},distinguishable_colors(3));

nexttile(3,[1 1])
showdat={mean(Prespike_ampSS(end-1:end,:),1,'omitnan')',mean(Prespike_ampCS(end-1:end,:),1,'omitnan')',mean(Prespike_ampBS(end-1:end,:),1,'omitnan')'};
plot_errorbar3([1:3],showdat,'ranksum',[-4 5],'Pre-spike AUC, distal dendrite',{'SS','CS','BS'},distinguishable_colors(3));

nexttile(4,[1 1])
showdat={mean(Prespike_ampSS,1,'omitnan')',mean(Prespike_ampCS,1,'omitnan')',mean(Prespike_ampBS,1,'omitnan')'};
plot_errorbar3([1:3],showdat,'ranksum',[-4 5],'Pre-spike AUC, all ROIs',{'SS','CS','BS'},distinguishable_colors(3));

figure(8); clf;
tiledlayout(2,2);
silentLap=unique(Result.VR(8,Result.Blue==0));
Unstimperiod=ismember(Result.VR(8,:),silentLap);
all_spike=find(max(SpClass([1 2 4],Unstimperiod),[],1));
Subthreshold=get_subthreshold(NormalizedTrace,all_spike,7,17);
BAsum=mean(Subthreshold(rois{1},Unstimperiod),1,'omitnan')+mean(Subthreshold(rois{2},Unstimperiod),1,'omitnan');
BAdiff=mean(Subthreshold(rois{1},Unstimperiod),1,'omitnan')-mean(Subthreshold(rois{2},Unstimperiod),1,'omitnan');

nexttile([1 1])
show_dat={(tovec(BAsum(all_spike'+[-5:-1]))),(tovec(BAsum(setdiff([1:length(Unstimperiod)],all_spike'+[-50:50]))))};
plot_errorbar3([1:2],show_dat,'ranksum',[-15 25],'Basal + Apical',{'Pre-spike','Silent period'},distinguishable_colors(3))
hold all
boxplot(cell2mat(show_dat'),[ones(1,length(show_dat{1})) ones(1,length(show_dat{2}))*2])
set(gca,'XTick',[1:2],'XTickLabel',{'Pre-spike','Silent period'})

nexttile([1 1])
show_dat={(tovec(BAdiff(all_spike'+[-5:-1]))),(tovec(BAdiff(setdiff([1:size(Result.traces,2)],all_spike'+[-50:50]))))};
plot_errorbar3([1:2],show_dat,'ranksum',[-0.5 15],' Basal - Apical ',{'Pre-spike','Silent period'},distinguishable_colors(3))
hold all
boxplot(cell2mat(show_dat'),[ones(1,length(show_dat{1})) ones(1,length(show_dat{2}))*2])
set(gca,'XTick',[1:2],'XTickLabel',{'Pre-spike','Silent period'})

ax1=nexttile([1 1]);
scatter_heatmap((tovec(BAdiff(all_spike'+[-5:-1]))),(tovec(BAsum(all_spike'+[-5:-1]))),50)
colormap(turbo)
xlabel('Basal - Apical'); ylabel('Basal + Apical')
ax2=nexttile([1 1]);
scatter_heatmap((tovec(BAdiff(setdiff([1:size(Result.traces,2)],all_spike'+[-50:50])))),(tovec(BAsum(setdiff([1:size(Result.traces,2)],all_spike'+[-50:50])))),50)
colormap(turbo)
xlabel('Basal - Apical'); ylabel('Basal + Apical')

linkaxes([ax1 ax2],'xy')

%% Compare SS and CS versus silent period
t_bin=[40:40:200];
figure(14); clf;
tiledlayout(2,length(t_bin))
loc=1; ax1=[]; cax=[-0.3 4];
for silent_time=t_bin
nTau_presilent={[-silent_time:40],[-silent_time:100],[-30:20]}; %SS, CS, dSP

% Isolated Somatic spike
SSpresilent_s=[];
SS_all=find(Result.SpClass(1,:));
for s=SS_all
    isnearby=sum(ismember(s+nTau_presilent{1},SS_all))>1;
    isnearbyCS=sum(ismember(s+nTau_presilent{1},find(Result.CStrace)))>1;
    if ~isnearby & ~isnan(s) & ~isnearbyCS
        SSpresilent_s=[SSpresilent_s s];
    end
end
sp_na=sum((SSpresilent_s'+nTau_presilent{1})<0 | (SSpresilent_s'+nTau_presilent{1})>size(Result.traces,2),2)==0;
SSpresilent_s=SSpresilent_s(sp_na);

% Isolated Complex spike
CSpresilent_s=[];
C_list=find(Result.SpClass(2,:));
CS_label=max(Result.spike,[],1).*bwlabel(Result.CStrace);

for g=1:length(C_list)
    s=C_list(g);
    s_tmp=Result.spike(1,:);
    s_tmp(find(CS_label==g))=0;
    isnearbyCS=max(bwlabel(Result.CStrace(s+nTau_presilent{2})))>1;
    isnearbyS=sum(ismember(s+nTau_presilent{2},find(s_tmp)))>0;
    if ~isnearbyCS & ~isnearbyS
        CSpresilent_s=[CSpresilent_s s];
    end
end
sp_na=sum((CSpresilent_s'+nTau_presilent{2})<0 | (CSpresilent_s'+nTau_presilent{2})>size(Result.traces,2),2)==0;
CSpresilent_s=CSpresilent_s(sp_na);

STA_CSmat_silent=squeeze(mean(reshape(NormalizedTrace_dirt(Result.dist_order(noi_dist),CSpresilent_s'+nTau_presilent{2}),sum(noi_dist),[],length(nTau_presilent{2})),2,'omitnan'));
STA_SSmat_silent=squeeze(mean(reshape(NormalizedTrace_dirt(Result.dist_order(noi_dist),SSpresilent_s'+nTau_presilent{1}),sum(noi_dist),[],length(nTau_presilent{1})),2,'omitnan'));
base_time=[1:silent_time];
STA_CSmat_silent=[STA_CSmat_silent-median(STA_CSmat_silent(:,base_time),2)];
STA_SSmat_silent=[STA_SSmat_silent-median(STA_SSmat_silent(:,base_time),2)];

ax1=[ax1 nexttile(loc,[1 1])];
cax_ss=[prctile(STA_SSmat_silent(:),10)  prctile(STA_SSmat_silent(:),95)];
imagesc(nTau_presilent{1},[1:sum(noi_dist)],STA_SSmat_silent,cax)
title(['SS with >' num2str(silent_time) 'ms silent period before' newline 'N = ' num2str(length(SSpresilent_s))])

ax1=[ax1 nexttile(loc+length(t_bin),[1 1])];
cax_cs=[prctile(STA_CSmat_silent(:),10)  prctile(STA_CSmat_silent(:),85)];
imagesc(nTau_presilent{2},[1:sum(noi_dist)],STA_CSmat_silent,cax)
title(['CS with >' num2str(silent_time) 'ms silent period before' newline 'N = ' num2str(length(CSpresilent_s))])
loc=loc+1;
end
linkaxes(ax1,'x')
colormap(turbo)
%% compare SS waveform upon blue

nROI=size(Result.normTraces,1);
spclass_ind=1;
%load aligned movie somatic spike
% alignmovlist=dir(fullfile(fpath{f},[alignedMovFN{spclass_ind} '*.tiff']));
% AlignMov=[];
% for l=1:length(alignmovlist)
%     l
%     AlignMov=cat(3,AlignMov,readtiff(fullfile(fpath{f},alignmovlist(l).name)));
% end
% sz_align=size(AlignMov);
% AlignMov=double(reshape(AlignMov,sz_align(1),sz_align(2),length(nTau{spclass_ind}),[]));
% AlignMov=AlignMov-mean(AlignMov(:,:,1:10,:),3);
% AlignMov=AlignMov(:,:,:,ismember(Result.StackedSpike{spclass_ind}(2,:),bAP_s));
% bAP_s(find(ismember(bAP_s,Result.StackedSpike{spclass_ind}(2,:))==0))=[];

bAP_s_on=find(Result.Blue(bAP_s)>0);
bAP_s_off=find(Result.Blue(bAP_s)==0);
% STAmov_on=mean(AlignMov(:,:,:,bAP_s_on),4);
% STAmov_off=mean(AlignMov(:,:,:,bAP_s_off),4);

figure(7); clf; cmap=distinguishable_colors(6); ax1=[]; tiledlayout(2,3)
cax=[-0.3 1];
ax1=[ax1 nexttile([1 1])];
STAss_Blueon=squeeze(mean(STA_SSmat(dist_order(noi_dist),bAP_s_on,:)-median(STA_SSmat(dist_order(noi_dist),bAP_s_on,1:-nTau{1}(1)),3,'omitnan'),2,'omitnan'));
STAss_Blueoff=squeeze(mean(STA_SSmat(dist_order(noi_dist),bAP_s_off,:)-median(STA_SSmat(dist_order(noi_dist),bAP_s_off,1:-nTau{1}(1)),3,'omitnan'),2,'omitnan'));
% STAss_Blueon=squeeze(mean(STA_SSmat(dist_order(noi_dist),bAP_s_on,:),2,'omitnan'));
% STAss_Blueoff=squeeze(mean(STA_SSmat(dist_order(noi_dist),bAP_s_off,:),2,'omitnan'));
%STAss_Blueon=STAss_Blueon-(mean(STAss_Blueon(:,-nTau{1}(1)+[-20:-5]),2)-mean(STAss_Blueoff(:,-nTau{1}(1)+[-20:-5]),2));
%STAss_Blueon=STAss_Blueon-(max(STAss_Blueon,[],2)-max(STAss_Blueoff,[],2));

imagesc(nTau{1},[1:length(noi)],STAss_Blueon,cax)
title('SS during blue on'); colormap(turbo);

ax1=[ax1 nexttile([1 1])];
imagesc(nTau{1},[1:length(noi)],STAss_Blueoff,cax)
title('SS during blue off'); colormap(turbo);

ax1=[ax1 nexttile([1 1])];
imagesc(nTau{1},[1:length(noi)],movmean(STAss_Blueon-STAss_Blueoff,1,2))
title('Blue On - Blue Off'); colormap(turbo);

ax1=[ax1 nexttile([1 1])];
l=plot(nTau{1},STAss_Blueon); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
ylim([-0.5 8])

ax1=[ax1 nexttile([1 1])];
l=plot(nTau{1},STAss_Blueoff); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
ylim([-0.5 8])

ax1=[ax1 nexttile([1 1])];
l=plot(nTau{1},movmean(STAss_Blueon-STAss_Blueoff,1,2)); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
linkaxes([ax1],'x')
xlim([-100 40])

%% compare CS waveform upon blue

nROI=size(Result.normTraces,1);

CS_s_on=find(Result.Blue(CS_s)>0);
CS_s_off=find(Result.Blue(CS_s)==0);

figure(8); clf; cmap=distinguishable_colors(6); ax1=[]; tiledlayout(2,3)
cax=[-0.5 1];
ax1=[ax1 nexttile([1 1])]; som_roi=find(dist_order(noi_dist)==1);
STAcs_Blueon=squeeze(mean(STA_CSmat(dist_order(noi_dist),CS_s_on,:)-median(STA_CSmat(dist_order(noi_dist),CS_s_on,1:-nTau{2}(1)),3,'omitnan'),2,'omitnan'));
STAcs_Blueoff=squeeze(mean(STA_CSmat(dist_order(noi_dist),CS_s_off,:)-median(STA_CSmat(dist_order(noi_dist),CS_s_off,1:-nTau{2}(1)),3,'omitnan'),2,'omitnan'));

%STAcs_Blueon=STAcs_Blueon/max(STAcs_Blueon(som_roi,-nTau{2}(1)+[0:2]),[],2);
%STAcs_Blueoff=STAcs_Blueoff/max(STAcs_Blueoff(som_roi,-nTau{2}(1)+[0:2]),[],2);
% STAcs_Blueon=squeeze(mean(STA_CSmat(dist_order(noi_dist),CS_s_on,:),2,'omitnan'));
% STAcs_Blueoff=squeeze(mean(STA_CSmat(dist_order(noi_dist),CS_s_off,:),2,'omitnan'));
%STAcs_Blueon=STAcs_Blueon-(mean(STAcs_Blueon(:,-nTau{2}(1)+[3:5]),2)-mean(STAcs_Blueoff(:,-nTau{2}(1)+[-50:-40]),2));
%STAcs_Blueon=STAcs_Blueon-(max(STAcs_Blueon,[],2)-max(STAcs_Blueoff,[],2));

imagesc(nTau{2},[1:length(noi)],STAcs_Blueon,cax)
title(['CS during blue on, N = ' num2str(length(CS_s_on))]); colormap(turbo);

ax1=[ax1 nexttile([1 1])];
imagesc(nTau{2},[1:length(noi)],STAcs_Blueoff,cax)
title(['CS during blue off, N = ' num2str(length(CS_s_off))]); colormap(turbo);

ax1=[ax1 nexttile([1 1])];
imagesc(nTau{2},[1:length(noi)],movmean(STAcs_Blueon-STAcs_Blueoff,3,2))
title('Blue On - Blue Off'); colormap(turbo);

ax1=[ax1 nexttile([1 1])];
l=plot(nTau{2},STAcs_Blueon); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
ylim([-1 7])

ax1=[ax1 nexttile([1 1])];
l=plot(nTau{2},STAcs_Blueoff); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
ylim([-1 7])

ax1=[ax1 nexttile([1 1])];
l=plot(nTau{2},movmean(STAcs_Blueon-STAcs_Blueoff,3,2)); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
linkaxes([ax1],'x')
xlim([-70 70])

%% compare BS waveform upon blue

nROI=size(Result.normTraces,1);

BS_s_on=find(Result.Blue(BS_s)>0);
BS_s_off=find(Result.Blue(BS_s)==0);

figure(9); clf; cmap=distinguishable_colors(6); ax1=[]; tiledlayout(2,3)
cax=[-0.5 2];
ax1=[ax1 nexttile([1 1])]; som_roi=find(dist_order(noi_dist)==1);
STAbs_Blueon=squeeze(mean(STA_BSmat(dist_order(noi_dist),BS_s_on,:)-median(STA_BSmat(dist_order(noi_dist),BS_s_on,1:-nTau{4}(1)),3,'omitnan'),2,'omitnan'));
STAbs_Blueoff=squeeze(mean(STA_BSmat(dist_order(noi_dist),BS_s_off,:)-median(STA_BSmat(dist_order(noi_dist),BS_s_off,1:-nTau{4}(1)),3,'omitnan'),2,'omitnan'));

%STAcs_Blueon=STAcs_Blueon/max(STAcs_Blueon(som_roi,-nTau{2}(1)+[0:2]),[],2);
%STAcs_Blueoff=STAcs_Blueoff/max(STAcs_Blueoff(som_roi,-nTau{2}(1)+[0:2]),[],2);
% STAcs_Blueon=squeeze(mean(STA_CSmat(dist_order(noi_dist),CS_s_on,:),2,'omitnan'));
% STAcs_Blueoff=squeeze(mean(STA_CSmat(dist_order(noi_dist),CS_s_off,:),2,'omitnan'));
%STAcs_Blueon=STAcs_Blueon-(mean(STAcs_Blueon(:,-nTau{2}(1)+[3:5]),2)-mean(STAcs_Blueoff(:,-nTau{2}(1)+[-50:-40]),2));
%STAcs_Blueon=STAcs_Blueon-(max(STAcs_Blueon,[],2)-max(STAcs_Blueoff,[],2));

imagesc(nTau{4},[1:length(noi)],STAbs_Blueon,cax)
title(['BS during blue on, N = ' num2str(length(BS_s_on))]); colormap(turbo);

ax1=[ax1 nexttile([1 1])];
imagesc(nTau{4},[1:length(noi)],STAbs_Blueoff,cax)
title(['BS during blue off, N = ' num2str(length(BS_s_off))]); colormap(turbo);

ax1=[ax1 nexttile([1 1])];
imagesc(nTau{4},[1:length(noi)],movmean(STAbs_Blueon-STAbs_Blueoff,3,2))
title('Blue On - Blue Off'); colormap(turbo);

ax1=[ax1 nexttile([1 1])];
l=plot(nTau{4},STAbs_Blueon); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
ylim([-1 7])

ax1=[ax1 nexttile([1 1])];
l=plot(nTau{4},STAbs_Blueoff); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
ylim([-1 7])

ax1=[ax1 nexttile([1 1])];
l=plot(nTau{4},movmean(STAbs_Blueon-STAbs_Blueoff,3,2)); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
linkaxes([ax1],'x')
xlim([-70 70])

%% Temporal correlation
rois={basal_ROI{f},apical_ROI{f},ref_ROI{f}}; autocorrBA=[]; crosscorrBA_Blue=[];
nTau_silent=[-5:100];
all_spike=find(max(Result.SpClass([1 2 3],:),[],1));
sp_na=sum((all_spike'+nTau_silent)<0 | (all_spike'+nTau_silent)>size(Result.traces,2),2)==0;
all_spike(~sp_na)=[];

BlueOnsetTime=find(Result.Blue>0,1);
SilentTrace=NormalizedTrace_dirt;
%SilentTrace(:,all_spike'+nTau_silent)=0; % Peri-spike to zero
SilentTrace(:,setdiff([1:size(Result.traces,2)],all_spike'+[-30:100]))=0; %Remain Peri-spikes
SilentTrace(:,Result.motionReject)=0;
SilentTrace(isnan(SilentTrace))=0;
SilentTrace=movmean(SilentTrace,5,2,'omitnan');

[SilentTrace_BA]=[mean(SilentTrace(rois{1},:),1); mean(SilentTrace(rois{2},:),1); mean(SilentTrace(rois{3},:),1)];
[autocorrBA(:,1)  lagTau]=xcorr(SilentTrace_BA(1,:),5000);
[autocorrBA(:,2)  lagTau]=xcorr(SilentTrace_BA(2,:),5000);
[autocorrBA(:,3)  lagTau]=xcorr(SilentTrace_BA(3,:),5000);
[crosscorrBA  lagTau]=xcorr(SilentTrace_BA(1,:),SilentTrace_BA(2,:),5000);

for r=1:3
    [autocorrBA_Blue{1}(:,r)  lagTau]=xcorr(SilentTrace_BA(r,1:BlueOnsetTime-1),5000); %{1} for Blue off, {2} for Blue On
    [autocorrBA_Blue{2}(:,r)  lagTau]=xcorr(SilentTrace_BA(r,BlueOnsetTime:end),5000);
end
[crosscorrBA_Blue{1}  lagTau]=xcorr(SilentTrace_BA(1,1:BlueOnsetTime-1),SilentTrace_BA(2,1:BlueOnsetTime-1),5000);
[crosscorrBA_Blue{2}  lagTau]=xcorr(SilentTrace_BA(1,BlueOnsetTime:end),SilentTrace_BA(2,BlueOnsetTime:end),5000);

figure(5); clf; tiledlayout(1,3); ax1=[];
ax1=[ax1 nexttile([1 1])];
plot(lagTau,crosscorrBA,'k'); hold all
plot(lagTau,autocorrBA); hold all
legend({'Cross Basal-Apical','Auto Basal','Auto Apical','Auto Soma'})
xlabel('Time (ms)');

ax1=[ax1 nexttile([1 1])];
plot(lagTau,crosscorrBA_Blue{1},'k'); hold all
plot(lagTau,autocorrBA_Blue{1}); hold all
legend({'Cross Basal-Apical','Auto Basal','Auto Apical','Auto Soma'})
title('Before Stimulation (start - Lap 12)'); xlabel('Time (ms)');

ax1=[ax1 nexttile([1 1])];
plot(lagTau,crosscorrBA_Blue{2},'k'); hold all
plot(lagTau,autocorrBA_Blue{2}); hold all
legend({'Cross Basal-Apical','Auto Basal','Auto Apical','Auto Soma'})
title('After Stimulation onset (Lap 13 - end)'); xlabel('Time (ms)');
linkaxes(ax1,'x')

%% moving correlation
window_sz=200; moving_bin=100; corrMat_moving=[]; corrAmp_moving=[]; nTau_silent=[-2:100];
t_slide=[1:moving_bin:size(Result.traces,2)];
for t=1:length(t_slide)-1
    t
    t_corr=[t_slide(t):t_slide(t)+window_sz];
    if t_slide(t)+window_sz<EndFrame(f)
    [corrMat_moving(t,:,1), t_lag]=xcorr(SilentTrace_BA(1,t_corr),SilentTrace_BA(2,t_corr),window_sz);
    T_center(t)=mean(t_corr);
    end
end
corrAmp_moving_int=interp1(T_center,corrMat_moving,[1:EndFrame(f)])';
all_spike=find(max(SpClass([1 2 4],:),[],1));
%corrAmp_moving_int(:,all_spike'+nTau_silent)=NaN;
corrAmp_moving_int(:,Result.motionReject)=NaN;
%corrAmp_moving_int(:,sum(Result.dirtTrace,1)>0)=NaN;

Corr_trace=mean(corrAmp_moving_int(window_sz+[-5:5],:),1,'omitnan');
Corr_Frnttrace=mean(corrAmp_moving_int(window_sz+[-120:-90],:),1,'omitnan');
Corr_Backtrace=mean(corrAmp_moving_int(window_sz+[90:120],:),1,'omitnan');

figure(23); clf;
tiledlayout(1,3);
nexttile([1 1])
show_dat={tovec(Corr_trace(all_spike'+[-5:-1])),tovec(Corr_trace(setdiff([1:length(Corr_trace)],all_spike'+[-50:50])))};
plot_errorbar3([1:2],show_dat,'ranksum',[-200 2000],'Correlation coefficient',{'Pre-spike','Silent period'},distinguishable_colors(3));

nexttile([1 1])
show_dat={abs(tovec(Corr_Frnttrace(all_spike'+[-5:-1]))),abs(tovec(Corr_Frnttrace(setdiff([1:length(Corr_trace)],all_spike'+[-50:50]))))};
plot_errorbar3([1:2],show_dat,'ranksum',[-200 1200],'Correlation coefficient',{'Pre-spike','Silent period'},distinguishable_colors(3))

nexttile([1 1])
show_dat={abs(tovec(Corr_Backtrace(all_spike'+[-5:-1]))),abs(tovec(Corr_Backtrace(setdiff([1:length(Corr_trace)],all_spike'+[-50:50]))))};
plot_errorbar3([1:2],show_dat,'ranksum',[-200 1200],'Correlation coefficient',{'Pre-spike','Silent period'},distinguishable_colors(3))


CorrAve_SSmat=reshape(Corr_trace(1,bAP_s(find(Result.Blue(bAP_s)==0))'+nTau{1}),1,[],length(nTau{1}));
CorrAve_CSmat=reshape(Corr_trace(1,CS_s(find(Result.Blue(CS_s)==0))'+nTau{2}),1,[],length(nTau{2}));
CorrAve_dSPmat=reshape(Corr_trace(1,dSP_s(find(Result.Blue(dSP_s)==0))'+nTau{3}),1,[],length(nTau{3}));
CorrAve_BSmat=reshape(Corr_trace(1,BS_s(find(Result.Blue(BS_s)==0))'+nTau{4}),1,[],length(nTau{4}));

CorrFrnt_SSmat=reshape(Corr_Frnttrace(1,bAP_s(find(Result.Blue(bAP_s)==0))'+nTau{1}),1,[],length(nTau{1}));
CorrFrnt_CSmat=reshape(Corr_Frnttrace(1,CS_s(find(Result.Blue(CS_s)==0))'+nTau{2}),1,[],length(nTau{2}));
CorrFrnt_dSPmat=reshape(Corr_Frnttrace(1,dSP_s(find(Result.Blue(dSP_s)==0))'+nTau{3}),1,[],length(nTau{3}));
CorrFrnt_BSmat=reshape(Corr_Frnttrace(1,BS_s(find(Result.Blue(BS_s)==0))'+nTau{4}),1,[],length(nTau{4}));

CorrBack_SSmat=reshape(Corr_Backtrace(1,bAP_s(find(Result.Blue(bAP_s)==0))'+nTau{1}),1,[],length(nTau{1}));
CorrBack_CSmat=reshape(Corr_Backtrace(1,CS_s(find(Result.Blue(CS_s)==0))'+nTau{2}),1,[],length(nTau{2}));
CorrBack_dSPmat=reshape(Corr_Backtrace(1,dSP_s(find(Result.Blue(dSP_s)==0))'+nTau{3}),1,[],length(nTau{3}));
CorrBack_BSmat=reshape(Corr_Backtrace(1,BS_s(find(Result.Blue(BS_s)==0))'+nTau{4}),1,[],length(nTau{4}));

Lap_corr = PlaceTrigger_average(corrAmp_moving_int,300,Result.VR,0.002,115); %total trace

figure(8); clf;
tiledlayout(4,4);
ax1=nexttile([1 4]);
imagesc([1:EndFrame(f)],t_lag,corrAmp_moving_int);
ax2=nexttile([1 4]);
imagesc(SilentTrace(Result.dist_order(noi_dist),:),[-1 3])
linkaxes([ax1 ax2],'x')
set(gcf,'WindowScrollWheelFcn',@scrollWheel)

nexttile([1 1])
Lap_CorrCenter=ringmovMean(mean(Lap_corr(:,:,find(t_lag<4 & t_lag>-4)),3,'omitnan'),5);
imagesc(Lap_CorrCenter)
title('Center correlation')

nexttile([1 1])
plot(mean(Lap_CorrCenter(StimOff_Lap,:),1,'omitnan')); hold all
plot(mean(Lap_CorrCenter(StimOn_Lap,:),1,'omitnan'))
legend({'Stim. off Lap','Stim. on Lap'});
xlabel('VR position (bin)')
title('Center correlation')

nexttile([1 1])
Lap_CorrEdge=ringmovMean(mean(Lap_corr(:,:,find(t_lag<-90 & t_lag>-120)),3,'omitnan') + mean(Lap_corr(:,:,find(t_lag>90 & t_lag<120)),3,'omitnan'),5);
imagesc(Lap_CorrEdge)
title('Edge correlation')

nexttile([1 1])
plot(mean(Lap_CorrEdge(StimOff_Lap,:),1,'omitnan')); hold all
plot(mean(Lap_CorrEdge(StimOn_Lap,:),1,'omitnan'))
legend({'Stim. off Lap','Stim. on Lap'})
xlabel('VR position (bin)')
title('Edge correlation')

nexttile([1 1])
plot(nTau{1},squeeze(mean(CorrAve_SSmat,2,'omitnan'))); hold all
plot(nTau{2},squeeze(mean(CorrAve_CSmat,2,'omitnan')))
plot(nTau{3},squeeze(mean(CorrAve_dSPmat,2,'omitnan')))
plot(nTau{4},squeeze(mean(CorrAve_BSmat,2,'omitnan')))
legend({'SS','CS','dSP','BS'}); xlabel('Peri-Spike time (ms)');

nexttile([1 1])
plot(nTau{1},squeeze(mean(CorrFrnt_SSmat,2,'omitnan'))); hold all
plot(nTau{2},squeeze(mean(CorrFrnt_CSmat,2,'omitnan')))
plot(nTau{3},squeeze(mean(CorrFrnt_dSPmat,2,'omitnan')))
plot(nTau{4},squeeze(mean(CorrFrnt_BSmat,2,'omitnan')))
legend({'SS','CS','dSP','BS'}); xlabel('Peri-Spike time (ms)');

nexttile([1 1])
plot(nTau{1},squeeze(mean(CorrBack_SSmat,2,'omitnan'))); hold all
plot(nTau{2},squeeze(mean(CorrBack_CSmat,2,'omitnan')))
plot(nTau{3},squeeze(mean(CorrBack_dSPmat,2,'omitnan')))
plot(nTau{4},squeeze(mean(CorrBack_BSmat,2,'omitnan')))
legend({'SS','CS','dSP','BS'}); xlabel('Peri-Spike time (ms)');

nexttile([1 1])
Corr_trace(isnan(Corr_trace))=0; Corr_Edgetrace(isnan(Corr_Edgetrace))=0;
[auto_CorrTrace tl]=xcorr(Corr_trace-mean(Corr_trace),10000);
[auto_CorrEdgeTrace tl]=xcorr(Corr_Edgetrace-mean(Corr_Edgetrace),10000);

plot(tl,auto_CorrTrace); hold all
plot(tl,auto_CorrEdgeTrace);
xlabel('Lag time (ms)');
legend({'Center corr','Edge corr'})

%% Place fields
velocity_threshold=0.002;
Lap_FR= PlaceTrigger_average(Result.spike(1,:),300,Result.VR,velocity_threshold,115); %total trace
Lap_sub= PlaceTrigger_average(Subthreshold,300,Result.VR,velocity_threshold,115); %total trace
[~, Lap_occupancy]= PlaceTrigger_average(double(Result.VR(end,:)>velocity_threshold),300,Result.VR,velocity_threshold,115); %total trace    
Lap_diff= PlaceTrigger_average(BAdiff,300,Result.VR,velocity_threshold,115); %total trace    

pre_PC_Lap=[2:6]; post_PC_Lap=[7:11]; StimOnLap=[12:21]; StimOffLap=[setdiff([1:max(Result.VR(8,:))],StimOnLap)];
%dendriteROIs={[36 37],[33 34],[31 32],[29 30 35],[1],[15],[2 3],[18 19 20],[8 9],[14 13],[4 5 6 7],[27 26],[22 23 24]};
dendriteROIs={Result.dist_order(noi_dist)};
cROIs=zeros(1,size(Result.ftprnt,3));
for d=1:length(dendriteROIs); cROIs(dendriteROIs{d})=d; end
cmap=distinguishable_colors(length(dendriteROIs));

figure(31); clf;
tiledlayout(2,2);
nexttile([1 1]);
imagesc(repmat(Lap_FR*1000,1,2)); hold all
scatter(ones(1,length(StimOnLap))*152,StimOnLap,30,[1 0 1],'filled','marker','>')
set(gca,'XTick',[1 150 300 450 600],'XTickLabel',[-2 -1 0 1 2])
xlabel('Position (m)')
ylabel('VR trials'); xlim([150 450]);
c=colorbar; c.Label.String = 'Firing rate (Hz)'; c.Label.Rotation = -90;
legend('Stimulation laps')

nexttile([1 1]);
imagesc(repmat(ringmovMean(Lap_diff,7),1,2)); hold all
scatter(ones(1,length(StimOnLap))*152,StimOnLap,30,[1 0 1],'filled','marker','>')
set(gca,'XTick',[1 150 300 450 600],'XTickLabel',[-2 -1 0 1 2])
title('Basal - Apical')
xlabel('Position (m)')
ylabel('VR trials'); xlim([150 450]);
c=colorbar; c.Label.String = 'Subthreshold'; c.Label.Rotation = -90;

nexttile([1 1])
Lap_prePC=ringmovMean(squeeze(mean(Lap_sub(pre_PC_Lap,:,cell2mat(dendriteROIs)),1,'omitnan'))',5);
imagesc(repmat((Lap_prePC),1,2)); hold all
scatter(ones(1,length(cell2mat(dendriteROIs)))*(152),[1:length(cell2mat(dendriteROIs))],40,cmap(cROIs(cell2mat(dendriteROIs)),:),'filled','marker','o')
set(gca,'XTick',[1 150 300 450 600],'XTickLabel',[-2 -1 0 1 2])
xlabel('Position (m)')
ylabel('Basal to apical dendrites'); xlim([150 450]);
c=colorbar; c.Label.String = 'Mean subthrehold'; c.Label.Rotation = -90;
title('Before place field formation, lap 3')

nexttile([1 1])
Lap_postPC=ringmovMean(squeeze(mean(Lap_sub(post_PC_Lap,:,cell2mat(dendriteROIs)),1,'omitnan'))',5);
imagesc(repmat((Lap_postPC),1,2)); hold all
scatter(ones(1,length(cell2mat(dendriteROIs)))*(152),[1:length(cell2mat(dendriteROIs))],40,cmap(cROIs(cell2mat(dendriteROIs)),:),'filled','marker','o')
set(gca,'XTick',[1 150 300 450 600],'XTickLabel',[-2 -1 0 1 2])
title('After place field formation, lap 5-11')
xlabel('Position (m)')
ylabel('Basal to apical dendrites'); xlim([150 450]);
c=colorbar; c.Label.String = 'Mean subthrehold'; c.Label.Rotation = -90;

figure(32); clf;
ind=cell2mat(dendriteROIs);
for n=1:length(ind)
nexttile([1 1])
lap_corr=[];
for l1=1:size(Lap_sub,1)
for l2=1:size(Lap_sub,1)
    nonnanbin=find(~isnan(Lap_sub(l1,:,ind(n))') & ~isnan(Lap_sub(l2,:,ind(n))'));
    if isempty(nonnanbin)
lap_corr(l1,l2)=NaN;
    else
        lap_corr(l1,l2)=corr(Lap_sub(l1,nonnanbin,ind(n))',Lap_sub(l2,nonnanbin,ind(n))');
    end
end
end
imagesc(lap_corr(StimOffLap,StimOffLap),[-0.3 0.6]); axis equal tight; hold all
title([num2str(n)],'color',cmap(cROIs(ind(n)),:))
xlabel('VR trials')
ylabel('VR trials')
end
colormap(turbo)
c=colorbar; title(c,'Correlation coefficient')

%% Clustering pre-spike dynamics
NclusterSSCS=[];
NormalizedTrace_dirt_HP=NormalizedTrace-movprc(NormalizedTrace,400,30,2);
NormalizedTrace_dirt_HP(Result.dirtTrace>0)=NaN;
NormalizedTrace_dirt_HP(:,Result.motionReject)=NaN;

STA_CSmat=reshape(NormalizedTrace_dirt_HP(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));
STA_SSmat=reshape(NormalizedTrace_dirt_HP(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1}));
STA_dSPmat=reshape(NormalizedTrace_dirt_HP(:,dSP_s'+nTau{3}),nROI,[],length(nTau{3}));
STA_BSmat=reshape(NormalizedTrace_dirt_HP(:,BS_s'+nTau{4}),nROI,[],length(nTau{4}));

STA_SSmat_on=STA_SSmat(:,bAP_s_on,:); STA_SSmat_off=STA_SSmat(:,bAP_s_off,:);%-mean(STA_SSmat(:,bAP_s_off,1:40),3,'omitnan');
STA_CSmat_on=STA_CSmat(:,CS_s_on,:); STA_CSmat_off=STA_CSmat(:,CS_s_off,:);%-mean(STA_CSmat(:,CS_s_off,1:40),3,'omitnan');
STA_BSmat_on=STA_BSmat(:,BS_s_on,:); STA_BSmat_off=STA_BSmat(:,BS_s_off,:);%-mean(STA_BSmat(:,BS_s_off,1:40),3,'omitnan');

 % PrespikeMat_SS=[tovec(permute(movmean(STA_SSmat_off(Result.dist_order(noi_dist),:,-nTau{1}(1)-10:-nTau{1}-2),6,3,'omitnan'),[1 3 2]))];
 % PrespikeMat_CS=[tovec(permute(movmean(STA_CSmat_off(Result.dist_order(noi_dist),:,-nTau{2}(1)-10:-nTau{2}-2),6,3,'omitnan'),[1 3 2]))];
 % PrespikeMat_BS=[tovec(permute(movmean(STA_BSmat_off(Result.dist_order(noi_dist),:,-nTau{4}(1)-10:-nTau{4}-2),6,3,'omitnan'),[1 3 2]))];
PrespikeMat_SS=[tovec(permute(mean(STA_SSmat_off(Result.dist_order(noi_dist),:,-nTau{1}(1)-6:-nTau{1}-2),3,'omitnan'),[1 3 2]))];
PrespikeMat_CS=[tovec(permute(mean(STA_CSmat_off(Result.dist_order(noi_dist),:,-nTau{2}(1)-6:-nTau{2}-2),3,'omitnan'),[1 3 2]))];
PrespikeMat_BS=[tovec(permute(mean(STA_BSmat_off(Result.dist_order(noi_dist),:,-nTau{4}(1)-6:-nTau{4}-2),3,'omitnan'),[1 3 2]))];

Vpre=[PrespikeMat_SS PrespikeMat_CS PrespikeMat_BS]; Classpre=[ones(1,size(PrespikeMat_SS,2)) ones(1,size(PrespikeMat_CS,2))*2 ones(1,size(PrespikeMat_BS,2))*3];
rmv_ind=sum(isnan(Vpre),1)>0;
Vpre(:,rmv_ind)=[]; Classpre(:,rmv_ind)=[];

[V,~, VprePC]=get_eigvector(Vpre,10); 
VprePC=VprePC';
%VprePC=Vpre;
%VprePC=rescale2(VprePC,1)';

totalN=[size(PrespikeMat_SS,2) size(PrespikeMat_CS,2) size(PrespikeMat_BS,2)];

%Vpre=(Vpre-prctile(Vpre,10,2))./((prctile(Vpre,90,2)-prctile(Vpre,10,2)));
distances = pdist(VprePC', 'euclidean');  % Euclidean distance
Z = linkage(distances, 'ward');
%num_clusters = 15;
cuttoff_val=22;
cluster_indices = cluster(Z, 'cutoff', cuttoff_val, 'criterion', 'distance');
num_clusters=max(cluster_indices);

Vpre_reduce = zeros(size(Vpre, 1), num_clusters); cluster_weight=zeros(1,num_clusters);

for i = 1:num_clusters
    Vpre_reduce(:, i) = mean(Vpre(:, find(cluster_indices == i)),2);
    %Vpre_reduce(:, i) = Vpre(:, find(cluster_indices == i,1));
    cluster_weight(i)=length(find(cluster_indices == i));
    Ncluster_class(i,:)=[sum(Classpre(find(cluster_indices == i))==1) sum(Classpre(find(cluster_indices == i))==2) sum(Classpre(find(cluster_indices == i))==3)];
end
Ncluster_class_fraction=Ncluster_class./totalN*100;
[~,sort_cluster]=sort(cluster_weight,'descend');
Vpre_reduce=Vpre_reduce(:,sort_cluster);
% Alternatively, you can plot a dendrogram to visualize the clustering
%imagesc(rescale2(Vs_reduce(:,sort_cluster),1))
figure(18); clf; tiledlayout(1,3);
nexttile([1 1]); dendrogram(Z,size(Z,1),'ColorThreshold',cuttoff_val);
nexttile([1 1]); silhouette(VprePC',cluster_indices,'Euclidean');
nexttile([1 1]); 
clusterBinMat=get_indMat(find_index_bh(sort_cluster,cluster_indices)');
for c=1:3
    scatter3(VprePC(1,find(clusterBinMat(:,c))),VprePC(2,find(clusterBinMat(:,c))),VprePC(3,find(clusterBinMat(:,c))),'filled'); hold all
end
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');

figure(15); clf; 
tiledlayout(1,num_clusters);
for n=1:num_clusters %Silent during blue off
    nexttile([1 1])
    imshow2(reshape(Vpre_reduce(:,n),sum(noi_dist),[]),[])
    colormap('turbo')
    title(['Cluster' num2str(n) ' , N = ' num2str(cluster_weight(sort_cluster(n))) newline ...
          'SS/CS/BS (%):' num2str(Ncluster_class_fraction(sort_cluster(n),1),2) '/' num2str(Ncluster_class_fraction(sort_cluster(n),2),2) '/' num2str(Ncluster_class_fraction(sort_cluster(n),3),2)])
end

figure(16); clf; %tiledlayout(5,num_clusters/5);
for n=1:num_clusters %Silent during blue off
    nexttile([1 1])
    meanV=mean(reshape(Vpre_reduce(:,n),sum(noi_dist),[]),2,'omitnan');
    meanV=meanV-min(meanV);
    imshow2(max(double(Result.ftprnt(:,:,Result.dist_order(noi_dist))>0).*reshape(meanV,1,1,[]),[],3),[]); colorbar;
    colormap('turbo')
    title(['Cluster' num2str(n) ' , N = ' num2str(cluster_weight(sort_cluster(n))) newline ...
          'SS/CS/BS (%) : ' num2str(Ncluster_class_fraction(sort_cluster(n),1),2) '/' num2str(Ncluster_class_fraction(sort_cluster(n),2),2) '/' num2str(Ncluster_class_fraction(sort_cluster(n),3),2)])
end

figure(17); clf; tiledlayout(3,8); cax=[-0.5 5]; ax1=[];
for n=1:num_clusters
    nexttile([1 1])
    meanV=mean(reshape(Vpre_reduce(:,n),sum(noi_dist),[]),2,'omitnan');
    meanV=meanV-min(meanV);
    imshow2(max(double(Result.ftprnt(:,:,Result.dist_order(noi_dist))>0).*reshape(meanV,1,1,[]),[],3),[]); colorbar;
    colormap('turbo')
    title(['Cluster' num2str(n) ' , N = ' num2str(cluster_weight(sort_cluster(n))) newline ...
          'SS/CS/BS (N) : ' num2str(Ncluster_class(sort_cluster(n),1),3) '/' num2str(Ncluster_class(sort_cluster(n),2),3) '/' num2str(Ncluster_class(sort_cluster(n),3),3)])
SS_ind=find(cluster_indices == sort_cluster(n) & Classpre' == 1); % SS 
CS_ind=find(cluster_indices == sort_cluster(n) & Classpre' == 2)-max(find(Classpre == 1)); % CS 
BS_ind=find(cluster_indices == sort_cluster(n) & Classpre' == 3)-max(find(Classpre == 2)); % BS
ax1=[ax1 nexttile([1 1])];
    imagesc(nTau{1},[1:sum(noi_dist)],squeeze(mean(STA_SSmat_off(Result.dist_order(noi_dist),SS_ind,:),2,'omitnan')),cax)
ax1=[ax1 nexttile([1 1])];
    imagesc(nTau{2},[1:sum(noi_dist)],squeeze(mean(STA_CSmat_off(Result.dist_order(noi_dist),CS_ind,:),2,'omitnan')),cax)
    ax1=[ax1 nexttile([1 1])];
    imagesc(nTau{4},[1:sum(noi_dist)],squeeze(mean(STA_BSmat_off(Result.dist_order(noi_dist),BS_ind,:),2,'omitnan')),cax)
end
linkaxes(ax1,'xy')
xlim([-30 10])

%% Correlation between dendrites
PrespikeMat_SS_ch=[]; PrespikeMat_CS_ch=[]; PrespikeMat_BS_ch=[]; PreMat=[];
STA_SSmat_off_ch=[]; STA_CSmat_off_ch=[]; STA_BSmat_off_ch=[];
for ch=1:2
STA_SSmat_off_ch{ch}=STA_SSmat_ch{ch}(:,bAP_s_off,:);
STA_CSmat_off_ch{ch}=STA_CSmat_ch{ch}(:,CS_s_off,:);
STA_BSmat_off_ch{ch}=STA_BSmat_ch{ch}(:,BS_s_off,:);

PrespikeMat_SS_ch{ch}=[(permute(movmean(STA_SSmat_off_ch{ch}(Result.dist_order(noi_dist),:,1:-nTau{1}-2),3,3,'omitnan'),[1 3 2]))];
PrespikeMat_CS_ch{ch}=[(permute(movmean(STA_CSmat_off_ch{ch}(Result.dist_order(noi_dist),:,1:-nTau{2}-2),3,3,'omitnan'),[1 3 2]))];
PrespikeMat_BS_ch{ch}=[(permute(movmean(STA_BSmat_off_ch{ch}(Result.dist_order(noi_dist),:,1:-nTau{4}-2),3,3,'omitnan'),[1 3 2]))];
PreMat{ch}=cat(3,PrespikeMat_SS_ch{ch},PrespikeMat_CS_ch{ch},PrespikeMat_BS_ch{ch}); 
end

SilentTrace_ch=cellfun(@(x) x(Result.dist_order(noi_dist),:),NormalizedTrace_ch,'UniformOutput',false);
SilentTrace_ch{1}(:,all_spike'+nTau_silent)=NaN; SilentTrace_ch{2}(:,all_spike'+nTau_silent)=NaN;
nonNan_frame=find(sum(isnan(SilentTrace_ch{1}),1)==0);
CorrMat=[];
for s=1:size(PreMat{ch},3)
for d1=1:size(PreMat{ch},1)
    for d2=1:size(PreMat{ch},1)
CorrMat(d1,d2,s)=corr(squeeze(PreMat{1}(d1,:,s))',squeeze(PreMat{2}(d2,:,s))');
    end
end
end

CorrMat_Silent=[];
for d1=1:size(PreMat{ch},1)
    for d2=1:size(PreMat{ch},1)
CorrMat_Silent(d1,d2)=corr(squeeze(SilentTrace_ch{1}(d1,nonNan_frame))',squeeze(SilentTrace_ch{2}(d2,nonNan_frame))');
    end
end

interDendDist=[];
SkelDend = Skeletonize_dendrite(Result.ref_im,6,0.01,10);
for i=1:nROI
i
for j=1:nROI
[interDendDist(i,j), path]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
end
end
interDendDist=interDendDist(Result.dist_order(noi_dist),Result.dist_order(noi_dist));

% moving time
window_sz=1000; moving_bin=200; corrMat_moving=[]; corrAmp_moving=[]; nTau_silent=[-2:100]; clustersIndicies=[];
k = 4; % Define number of clusters (k)
CorrMat_silentMoving=[]; ARI=[]; T_center=[]; CorrMat_silentMoving_avg=[];
t_slide=[1:moving_bin:size(Result.traces,2)]; g=1;

Ct = squeeze(CorrMat_Silent);
Ct(logical(eye(size(Ct)))) = 1;
distMatrix = sqrt(2*(1 - Ct));
Z_ref = linkage(squareform(distMatrix), 'ward');
% distMatrix=pdist(SilentTrace_ch{1}(:,nonNan_frame),"euclidean");
% Z_ref = linkage((distMatrix), 'ward');
leafOrder = optimalleaforder(Z_ref,squareform(distMatrix));
Cluster_ref= switchlabel(cluster(Z_ref, 'maxclust', k));
[~, order_cluster]=sort(Cluster_ref,'ascend');

for t=1:length(t_slide)-1
    t
    t_corr=setdiff([t_slide(t):t_slide(t)+window_sz],find(sum(isnan(SilentTrace_ch{1}),1)>0)); %remove NaN frame
    if ~isempty(t_corr) & t_slide(t)+window_sz<EndFrame(f)

        for d1=1:size(PreMat{ch},1)
            for d2=1:size(PreMat{ch},1)
                CorrMat_silentMoving(d1,d2,g)=corr(squeeze(SilentTrace_ch{1}(d1,t_corr))',squeeze(SilentTrace_ch{2}(d2,t_corr))');
            end
        end
        T_center(g)=mean(t_corr);

       Ct = squeeze(CorrMat_silentMoving(:, :, g));
       Ct(logical(eye(size(Ct)))) = 1;
    distMatrix = sqrt(2*(1 - Ct)); % This is one possible way to define distance (1 - correlation)
    Z = linkage(squareform(distMatrix), 'ward');
    clustersIndicies(:,g) = cluster(Z, 'maxclust', k);
    clustersIndicies(:,g)=switchlabel(clustersIndicies(:,g));
    ARI(g)=jaccardSimilarity(Cluster_ref, clustersIndicies(:,g));
    g=g+1;

    end
end

for cl1=1:3
for cl2=1:3
    CorrMat_silentMoving_avg(cl1,cl2,:)=mean(CorrMat_silentMoving(Cluster_ref==cl1,Cluster_ref==cl2,:),[1 2]);
end
end

figure(22); clf; ax_corr=[];
tiledlayout(4,6)
nexttile([1 3])
show_footprnt_contour(Result.ftprnt(:,:,Result.dist_order(noi_dist)),Result.ref_im);
%colormap(ax_fprnt,'gray')

ax_fprnt=nexttile(4,[1 3]);
imagesc(max(double(Result.ftprnt(:,:,Result.dist_order(noi_dist))>0).*reshape(Cluster_ref,1,1,[]),[],3)); axis equal tight off
colormap(ax_fprnt,'sky')
title('Classified cluster')

ax_corr=nexttile([1 1]);
imagesc(mean(CorrMat,3,'omitnan'),[0 0.7]); axis tight equal
title('Correlation during Pre-spike time')
colormap(ax_corr,'turbo')

ax_corr=nexttile([1 1]);
imagesc(CorrMat_Silent,[0 0.7]); axis tight equal
title('Correlation during silent')
colormap(ax_corr,'turbo')

ax_corr=nexttile([1 1]);
imagesc(CorrMat_Silent(order_cluster,order_cluster),[0 0.7]); axis tight equal
title('Correlation during silent, reorder by cluster')
colormap(ax_corr,'turbo')

nexttile([1 2]);
dendrogram(Z_ref,0,'ColorThreshold',1.45); axis off

ax_corr=nexttile([1 6]);
imagesc(T_center,[1:sum(noi_dist)],clustersIndicies)
title('3 Cluster over silent time')
colormap(ax_corr,'sky')

nexttile([1 2])
for cl=1:3
plot(T_center,squeeze(CorrMat_silentMoving_avg(cl,cl,:))); hold all
end
xlabel('Time (ms)')
ylabel('Correlation coefficient')
legend({'Cluster 1','Cluster 2','Cluster 3'});

nexttile([1 2])
plot(T_center,squeeze(CorrMat_silentMoving_avg(1,2,:))); hold all
plot(T_center,squeeze(CorrMat_silentMoving_avg(1,3,:))); hold all
plot(T_center,squeeze(CorrMat_silentMoving_avg(2,3,:))); hold all
xlabel('Time (ms)')
ylabel('Correlation coefficient')
legend({'Cluster 1 & 2','Cluster 1 & 3','Cluster 2 & 3'});

nexttile([1 2])
for cl=1:k
ind=triu(CorrMat_Silent,1)>0 & ((Cluster_ref==cl)*(Cluster_ref==cl)');
scatter(tovec(interDendDist(ind))*1.17,tovec(CorrMat_Silent(ind)),'filled'); hold all
end
legend({'cluster 1','cluster 2','cluster 3'})
xlabel('Pairwise geodesic distance (\mum)')
ylabel('Correlation coefficient')
%% Clustering Pre-spike dynamics in movie



%% Blue light unmasked subthreshold place field
nBin=150; Lap_sub=[]; StimOn_Lap=[12:21];
all_spike=find(max(Result.SpClass([1 2],:),[],1));
Subthreshold=get_subthreshold(NormalizedTrace,all_spike,7,17);
Subthreshold_silent=Subthreshold;
sp_na=sum((all_spike'+[-8:20])<0 | (all_spike'+[-8:20])>size(Result.traces,2),2)==0;
all_spike=all_spike(sp_na);
Subthreshold_silent(:,all_spike'+[-8:20])=NaN;
Subthreshold_silent(:,Result.motionReject)=NaN;
Subthreshold_silent(:,Result.CStrace>0)=NaN;
Subthreshold_silent(Result.dirtTrace>0)=NaN;
bwBlue=bwlabel(Result.Blue>0);

for r=1:3 %basal, apical, soma
    Lap_sub{r}=PlaceTrigger_average(mean(Subthreshold_silent(rois{r},:),1,'omitnan'),nBin,Result.VR,0,115);
end
Lap_FR = PlaceTrigger_average(Result.spike(1,:),nBin,Result.VR,0,115); %total trace
Lap_CS = PlaceTrigger_average(Result.SpClass(2,:),nBin,Result.VR,0,115); %total trace
Sp_tr=Result.spike(1,:); Sp_tr(Result.Blue==0)=NaN; %during blue on
Lap_FR_BlueOn=PlaceTrigger_average(Sp_tr,nBin,Result.VR,0,115); %soma during blue on
Lap_blue=PlaceTrigger_average(Result.Blue,nBin,Result.VR,0,115); %Blue

figure(10); clf; tiledlayout(3,3); ax1=[]; ax3=[];
ax1=[ax1 nexttile([1 1])];
img=imagesc([1:nBin]*2/nBin,[1:size(Lap_FR,1)],ringmovMean(Lap_FR,5)); hold all
title('Firing rate'); xlabel('Position (m)'); ylabel('Laps');
ax1=[ax1 nexttile([1 1])];
img2=imagesc([1:nBin]*2/nBin,[1:size(Lap_FR,1)],ringmovMean(Lap_CS,5));
arrayfun(@(x) colormap(x,"turbo"),ax1);
title('Complex spike rate'); xlabel('Position (m)'); ylabel('Laps');
ax2=[nexttile([1 1])];
img3=imagesc([1:nBin]*2/nBin,[1:size(Lap_FR,1)],Lap_blue>0); hold all
arrayfun(@(x) colormap(x,"gray"),ax2)
title('Blue Stimulation'); xlabel('Position (m)'); ylabel('Laps');
for r=1:3
    ax3=[ax3 nexttile([1 1])];
    imagesc(ringmovMean(Lap_sub{r},5)); hold all
    plot(ones(1,length(StimOn_Lap))*3,StimOn_Lap,'marker','>','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[0 0 0])
    title([title_str{r} ' subthreshold'])
end
legend({'Stimulation Lap'},'Location','northwest')
arrayfun(@(x) colormap(x,gen_colormap([0 0.5 1; 1 1 1; 1 0 0])),ax3)

nexttile([1 1]);
Lap_sub_mean=cell2mat(cellfun(@(x) mean(ringmovMean(x(3:12,:),5),1,'omitnan'),Lap_sub,'UniformOutput',false)');
plot([1:nBin]*2/nBin,Lap_sub_mean')
legend(title_str);
xlabel('VR Position (m)'); ylabel('Subthreshold (A.U.)');
title('Mean subthreshold during first 10 laps')

nexttile([1 1]);
plot([1:nBin]*2/nBin,mean(ringmovMean(Lap_FR_BlueOn,5),1,'omitnan')*1000);
xlabel('VR Position (m)'); ylabel('Firing rate (Hz)');
title('Mean firing rate during stimulation')

nexttile([1 1]);
scatter(mean(ringmovMean(Lap_FR_BlueOn,5),1,'omitnan')*1000,Lap_sub_mean,'filled')
[corrval pval]=corr(mean(ringmovMean(Lap_FR_BlueOn,5),1,'omitnan')',Lap_sub_mean');
legend(strcat(title_str', {', corr : '}, num2str(corrval',2)),'Location','southeast');
box on; ylabel('Subthreshold (A.U.)'); xlabel('Firing rate (Hz)');
title('Correlation btw Subth. and Firing rate')

%% Theta osillation, phase, amplitude
theta_trace=[]; phase_trace=[]; phaseMag_trace=[];
Hilbert_SSmat=[]; Hilbert_CSmat=[];
theta_pass_frq=[5 12]; 
freq_lowhigh=theta_pass_frq/(1000/2);
[b4, a4] = butter(4, freq_lowhigh, 'bandpass');

for n=1:nROI
    subthreshold_rmvNaN=Subthreshold(n,:);
    subthreshold_rmvNaN(isnan(subthreshold_rmvNaN))=median(subthreshold_rmvNaN,'omitnan');
    theta_trace(n,:)=filtfilt(b4, a4, subthreshold_rmvNaN);
    phase_trace(n,:) = angle(hilbert(theta_trace(n,:)));
    phaseMag_trace(n,:)=abs(hilbert(theta_trace(n,:)));
end

CS_s_new=find(Result.SpClass(2,:)); CS_s_new(Result.Blue(CS_s_new)>0 | sum(isnan(reshape(NormalizedTrace_dirt(:,CS_s_new'+nTau{2}),nROI,[],length(nTau{2}))),[1 3])>0)=[]; % spikes during blue off
SS_s_new=bAP_s; SS_s_new(Result.Blue(SS_s_new)>0)=[]; 
dSP_s_new=dSP_s; dSP_s_new(Result.Blue(dSP_s_new)>0)=[]; 

theta_trace_dirt=theta_trace;
theta_trace_dirt(Result.dirtTrace>0)=NaN;
theta_trace_dirt(:,Result.motionReject)=NaN;

phase_trace_dirt=phase_trace;
phase_trace_dirt(Result.dirtTrace>0)=NaN;
phase_trace_dirt(:,Result.motionReject)=NaN;

phaseMag_trace_dirt=phaseMag_trace;
phaseMag_trace_dirt(Result.dirtTrace>0)=NaN;
phaseMag_trace_dirt(:,Result.motionReject)=NaN;

Stheta_CSmat=reshape(theta_trace_dirt(:,CS_s_new'+nTau{2}),nROI,[],length(nTau{2}));
Stheta_SSmat=reshape(theta_trace_dirt(:,SS_s_new'+nTau{1}),nROI,[],length(nTau{1}));
Stheta_dSPmat=reshape(theta_trace_dirt(:,dSP_s_new'+nTau{3}),nROI,[],length(nTau{3}));

SthetaPhase_CSmat=reshape(phase_trace_dirt(:,CS_s_new'+nTau{2}),nROI,[],length(nTau{2}));
SthetaPhase_SSmat=reshape(phase_trace_dirt(:,SS_s_new'+nTau{1}),nROI,[],length(nTau{1}));
SthetaPhase_dSPmat=reshape(phase_trace_dirt(:,dSP_s_new'+nTau{3}),nROI,[],length(nTau{3}));

SthetaMag_CSmat=reshape(phaseMag_trace_dirt(:,CS_s_new'+nTau{2}),nROI,[],length(nTau{2}));
SthetaMag_SSmat=reshape(phaseMag_trace_dirt(:,SS_s_new'+nTau{1}),nROI,[],length(nTau{1}));
SthetaMag_dSPmat=reshape(phaseMag_trace_dirt(:,dSP_s_new'+nTau{3}),nROI,[],length(nTau{3}));

figure(50); clf; cmap=distinguishable_colors(6); ax1=[];
noi=setdiff([1:nROI],[BadROI{f}]);
noi_dist=ismember(Result.dist_order,noi);
%rois={basal_ROI{f},apical_ROI{f}};
cax=[-1 1];%*0.01;
tiledlayout(3,6);
ax3=nexttile(1,[1 1]);
imagesc(nTau{1},[1:length(noi)],squeeze(mean(Stheta_SSmat(Result.dist_order(noi_dist),:,:),2,'omitnan')),cax)
title('Simple spike')
ax3=[ax3 nexttile(3,[1 1])];
imagesc(nTau{2},[1:length(noi)],squeeze(mean(Stheta_CSmat(Result.dist_order(noi_dist),:,:),2,'omitnan')),cax)
title('Complex spike')
ax3=[ax3 nexttile(5,[1 1])];
%imagesc(nTau{3},[1:length(noi)],squeeze(mean(Stheta_dSPmat(Result.dist_order(noi_dist),:,:),2,'omitnan')),cax)
%title('Dendritic spike')
%colormap(turbo); linkaxes(ax3,'xy'); xlim([-140 0])

ax1=[ax1 nexttile(2,[1 1])];
l=plot(nTau{1},squeeze(mean(Stheta_SSmat(Result.dist_order(noi_dist),:,:),2,'omitnan'))); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
xlabel('Peri-spike time (ms)'); ylabel('Theta filtered trace');
ax1=[ax1 nexttile(4,[1 1])];
l=plot(nTau{2},squeeze(mean(Stheta_CSmat(Result.dist_order(noi_dist),:,:),2,'omitnan'))); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
xlabel('Peri-spike time (ms)'); ylabel('Theta filtered trace');
ax1=[ax1 nexttile(6,[1 1])];
l=plot(nTau{3},squeeze(mean(Stheta_dSPmat(Result.dist_order(noi_dist),:,:),2,'omitnan'))); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
xlabel('Peri-spike time (ms)'); ylabel('Theta filtered trace');
linkaxes([ax1],'xy'); xlim([-140 0])

ax4=nexttile([2 2]); % Phase plot SS
l=polarplot(squeeze(mean_angle(SthetaPhase_SSmat(Result.dist_order(noi_dist),:,1:-nTau{1}(1)+1),2))',squeeze(mean(SthetaMag_SSmat(Result.dist_order(noi_dist),:,1:-nTau{1}(1)+1),2,'omitnan'))');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
hold all
polarplot(squeeze(mean_angle(SthetaPhase_SSmat(Result.dist_order(noi_dist),:,1),2))',squeeze(mean(SthetaMag_SSmat(Result.dist_order(noi_dist),:,1),2,'omitnan'))','k>');
polarplot(squeeze(mean_angle(SthetaPhase_SSmat(Result.dist_order(noi_dist),:,-nTau{1}(1)+1),2))',squeeze(mean(SthetaMag_SSmat(Result.dist_order(noi_dist),:,-nTau{1}(1)+1),2,'omitnan'))','ro');
rlim([0, 1.8]); title('SS')

ax4=nexttile([2 2]); % Phase plot CS
l=polarplot(squeeze(mean_angle(SthetaPhase_CSmat(Result.dist_order(noi_dist),:,1:-nTau{2}(1)+1),2))',squeeze(mean(SthetaMag_CSmat(Result.dist_order(noi_dist),:,1:-nTau{2}(1)+1),2,'omitnan'))');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
hold all
polarplot(squeeze(mean_angle(SthetaPhase_CSmat(Result.dist_order(noi_dist),:,1),2))',squeeze(mean(SthetaMag_CSmat(Result.dist_order(noi_dist),:,1),2,'omitnan'))','k>');
polarplot(squeeze(mean_angle(SthetaPhase_CSmat(Result.dist_order(noi_dist),:,-nTau{2}(1)+1),2))',squeeze(mean(SthetaMag_CSmat(Result.dist_order(noi_dist),:,-nTau{2}(1)+1),2,'omitnan'))','ro');
rlim([0, 1.8]); title('CS')

ax4=nexttile([2 2]); % Phase plot dSP
l=polarplot(squeeze(mean_angle(SthetaPhase_dSPmat(Result.dist_order(noi_dist),:,1:-nTau{3}(1)+1),2))',squeeze(mean(SthetaMag_dSPmat(Result.dist_order(noi_dist),:,1:-nTau{3}(1)+1),2,'omitnan'))');
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
hold all
polarplot(squeeze(mean_angle(SthetaPhase_dSPmat(Result.dist_order(noi_dist),:,1),2))',squeeze(mean(SthetaMag_dSPmat(Result.dist_order(noi_dist),:,1),2,'omitnan'))','k>');
polarplot(squeeze(mean_angle(SthetaPhase_dSPmat(Result.dist_order(noi_dist),:,-nTau{3}(1)+1),2))',squeeze(mean(SthetaMag_dSPmat(Result.dist_order(noi_dist),:,-nTau{3}(1)+1),2,'omitnan'))','ro');
rlim([0, 1.8]); title('dSP')

figure(52); clf; tiledlayout(3,3)
angle_bin=[-pi:0.3:pi];
ax5=nexttile([1 1]); % Phase histogram SS
histogram(squeeze(mean_angle(SthetaPhase_SSmat(rois{1},:,-nTau{1}(1)),1)),angle_bin,'Normalization','probability','FaceColor',[0 0 0]); hold all
histogram(squeeze(mean_angle(SthetaPhase_CSmat(rois{1},:,-nTau{2}(1)),1)),angle_bin,'Normalization','probability','FaceColor',[1 0 0]); hold all
histogram(squeeze(mean_angle(SthetaPhase_dSPmat(rois{1},:,-nTau{3}(1)),1)),angle_bin,'Normalization','probability','FaceColor',[0 0.5 1]); hold all
xlabel('Phase'); ylabel('Probability'); title('Basal phase');
legend({'SS','CS','dSp'})

ax5=nexttile([1 1]); % Phase histogram CS
histogram(squeeze(mean_angle(SthetaPhase_SSmat(rois{2},:,-nTau{1}(1)),1)),angle_bin,'Normalization','probability','FaceColor',[0 0 0]); hold all
histogram(squeeze(mean_angle(SthetaPhase_CSmat(rois{2},:,-nTau{2}(1)),1)),angle_bin,'Normalization','probability','FaceColor',[1 0 0]); hold all
histogram(squeeze(mean_angle(SthetaPhase_dSPmat(rois{2},:,-nTau{3}(1)),1)),angle_bin,'Normalization','probability','FaceColor',[0 0.5 1]); hold all
xlabel('Phase'); ylabel('Probability'); title('Apical phase');
legend({'SS','CS','dSp'})

ax5=nexttile([1 1]); % Phase histogram dSP
histogram(squeeze(mean_angle(SthetaPhase_SSmat(rois{3},:,-nTau{1}(1)),1)),angle_bin,'Normalization','probability','FaceColor',[0 0 0]); hold all
histogram(squeeze(mean_angle(SthetaPhase_CSmat(rois{3},:,-nTau{2}(1)),1)),angle_bin,'Normalization','probability','FaceColor',[1 0 0]); hold all
histogram(squeeze(mean_angle(SthetaPhase_dSPmat(rois{3},:,-nTau{3}(1)),1)),angle_bin,'Normalization','probability','FaceColor',[0 0.5 1]); hold all
xlabel('Phase'); ylabel('Probability'); title('Peri-soma phase');
legend({'SS','CS','dSp'}) 

ax1=[ax1 nexttile([1 1])];
cmap=[0 0 0;1 0 0;0 0.6 1]; % Position, Phase plot
for c=1:3
plot(repmat(Result.VR(5,find(Result.SpClass(c,:) & Result.Blue==0)),3,1),phase_trace(ref_ROI{f}(1),find(Result.SpClass(c,:) & Result.Blue==0))+[-pi 0 pi]','.','color',cmap(c,:),'MarkerSize',12); hold all
plot(repmat(Result.VR(5,find(Result.SpClass(c,:) & Result.Blue==0)),3,1)+115,phase_trace(ref_ROI{f}(1),find(Result.SpClass(c,:) & Result.Blue==0))+[-pi 0 pi]','.','color',cmap(c,:),'MarkerSize',12); hold all
end
xlabel('VR position'); ylabel('Phase')
ylim([-pi*3/4 pi*3/4]); xlim([70 150]);

ax1=[ax1 nexttile([1 1])];
errorbar_shade(nTau{1},squeeze(mean(mean(SthetaMag_SSmat(rois{1},:,:),1,'omitnan'),2,'omitnan')),squeeze(std(mean(SthetaMag_SSmat(rois{1},:,:),1,'omitnan'),0,2,'omitnan')),[0 0 0]); hold all
errorbar_shade(nTau{2},squeeze(mean(mean(SthetaMag_CSmat(rois{1},:,:),1,'omitnan'),2,'omitnan')),squeeze(std(mean(SthetaMag_CSmat(rois{1},:,:),1,'omitnan'),0,2,'omitnan')),[1 0 0]); hold all
xlabel('Peri-spike time (ms)'); ylabel('Theta magnitude'); title('Basal theta');
xlim([-140 0])

ax1=[ax1 nexttile([1 1])];
errorbar_shade(nTau{1},squeeze(mean(mean(SthetaMag_SSmat(rois{2},:,:),1,'omitnan'),2,'omitnan')),squeeze(std(mean(SthetaMag_SSmat(rois{2},:,:),1,'omitnan'),0,2,'omitnan')),[0 0 0]); hold all
errorbar_shade(nTau{2},squeeze(mean(mean(SthetaMag_CSmat(rois{2},:,:),1,'omitnan'),2,'omitnan')),squeeze(std(mean(SthetaMag_CSmat(rois{2},:,:),1,'omitnan'),0,2,'omitnan')),[1 0 0]); hold all
xlabel('Peri-spike time (ms)'); ylabel('Theta magnitude'); title('Apical theta');
xlim([-140 0])

for r=[1:3]
t=[-200:-30];    
nexttile([1 1])
    for c=1:3
sp=find(Result.SpClass(c,:) & Result.Blue==0);
sp(find(sp-120<0 | sp+50>size(Result.traces,2)))=[];
rng_theta=max(reshape(mean(theta_trace_dirt(rois{r},sp'+t),1,'omitnan'),[],length(t)),[],2)-min(reshape(mean(theta_trace_dirt(rois{r},sp'+t),1,'omitnan'),[],length(t)),[],2);
plot(rng_theta,mean_angle(phase_trace_dirt(rois{r},sp),1),'.','color',cmap(c,:),'MarkerSize',12); hold all
    end
    xlabel(['Theta swing before spike']);
    ylabel(['Phase at spike'])
    title([title_str{r}])
end

%% Upward going and downward going
binEdge=[-5:0.1:5]; up_thres=2; dw_thres=1.7; smooth_window=3000;
blueTargetF=mean(reshape(NormalizedTrace(:,bAP_ref'+[6:10]),size(NormalizedTrace,1),length(bAP_ref),[]),[2 3]);

bwBlue=bwlabel(Result.Blue>0);
all_spike=find(max(SpClass([1 2 4],:),[],1));
Subthreshold=get_subthreshold(NormalizedTrace,all_spike,7,17);
Subthreshold_silent=Subthreshold;
sp_na=sum((all_spike'+[-8:20])<0 | (all_spike'+[-8:20])>size(Result.traces,2),2)==0;
all_spike=all_spike(sp_na);
Subthreshold_silent(:,all_spike'+[-8:20])=NaN;
Subthreshold_silent(:,Result.motionReject)=NaN;
Subthreshold_silent(:,Result.CStrace>0)=NaN;
Subthreshold_silent(Result.dirtTrace>0)=NaN;
%SubSilent_hipass=movmean(Subthreshold_silent,10,2)-movprc(Subthreshold_silent,smooth_window,30,2);
SubSilent_hipass_up=movmean(Subthreshold_silent,10,2,'omitnan');%-movmedian(Subthreshold_silent,smooth_window,2,'omitnan');
SubSilent_hipass_down=blueTargetF-movmean(Subthreshold_silent,10,2,'omitnan');%-movmedian(Subthreshold_silent,smooth_window,2,'omitnan');

MedianF=zeros(2,max(bwBlue),3); 
MedianVRposition=zeros(2,max(bwBlue)); 

[upward_trace]=detect_transientSimple(SubSilent_hipass,[up_thres 1]);
[downward_trace]=detect_transientSimple(blueTargetF-SubSilent_hipass,[dw_thres 1]);

[~, ~, BlueOnErode]=get_blueoffTrace(double(Result.Blue>0),Result.Blue>0,100,100);
upward_trace(:,~BlueOnErode)=0;
[~, BlueErode]=get_blueoffTrace(double(Result.Blue>0),Result.Blue==0,100,500);
downward_trace(:,~BlueErode)=0;

figure(24); clf;
tiledlayout(4,1)
ax1=nexttile([1 1]);
imagesc(Subthreshold(cell2mat(dendriteROIs),:),[-1 4])
colormap(turbo)
ax2=nexttile([1 1]);
imagesc(upward_trace(cell2mat(dendriteROIs),:)>0)
ax3=nexttile([1 1]);
imagesc(downward_trace(cell2mat(dendriteROIs),:)>0)
ax4=nexttile([1 1]);
plot(Result.Blue)
linkaxes([ax1 ax2 ax3 ax4],'x')
set(gcf,'WindowScrollWheelFcn',@scrollWheel)




    for b=1:max(bwBlue)
        t_blue=[find(bwBlue==b,1)+100:max(find(bwBlue==b))];
        %dSub=mean(Subthreshold_silent(rois{r},t_blue),1)-movmedian(mean(Subthreshold_silent(rois{r},t_blue),1,'omitnan'),smooth_window,'omitnan');
        dSub=mean(Subthreshold_silent(rois{r},t_blue),1)-blueTargetF;
        dSub=movmean(dSub,10);
        MedianF(1,b,r)=median(mean(Subthreshold_silent(rois{r},t_blue),1,'omitnan'),'omitnan'); MedianVRposition(1,b)=median(Result.VR(5,t_blue));

        [Pk,pk_time]=findpeaks(dSub,'MinPeakDistance',10);
        [Pk_inv,pkinv_time]=findpeaks(-dSub,'MinPeakDistance',10);
        upward_pk{r,1}=[upward_pk{r,1} Pk];
        downward_pk{r,1}=[downward_pk{r,1} -Pk_inv];

        t_blueoff=t_blue(300:end-300)+2000;
        dSubOff=mean(Subthreshold_silent(rois{r},t_blueoff),1)-movmedian(mean(Subthreshold_silent(rois{r},t_blueoff),1,'omitnan'),smooth_window,'omitnan');
        dSubOff=movmean(dSubOff,5);
        MedianF(2,b,r)=median(mean(Subthreshold_silent(rois{r},t_blueoff),1,'omitnan'),'omitnan'); MedianVRposition(2,b)=median(Result.VR(5,t_blueoff));

        [Pk,pk_time]=findpeaks(dSubOff,'MinPeakDistance',10);
        [Pk_inv,pkinv_time]=findpeaks(-dSubOff,'MinPeakDistance',10);
        upward_pk{r,2}=[upward_pk{r,2} Pk];
        downward_pk{r,2}=[downward_pk{r,2} -Pk_inv];

        [~, up_traceBon]=detect_transient(dSub,[up_thres 1],zeros(1,length(t_blue)));
        [~, dw_traceBon]=detect_transient(-dSub,[dw_thres 1],zeros(1,length(t_blue)));
        [~, up_traceBoff]=detect_transient(dSubOff,[up_thres 1],zeros(1,length(t_blueoff)));
        [~, dw_traceBoff]=detect_transient(-dSubOff,[dw_thres 1],zeros(1,length(t_blueoff)));

        upward_trace(r,t_blue)=up_traceBon; upward_trace(r,t_blueoff)=up_traceBoff;
        downward_trace(r,t_blue)=dw_traceBon; downward_trace(r,t_blueoff)=dw_traceBoff;
    end

    % after blue
    t=[max(find(bwBlue>0))+1:size(Result.traces,2)];
    dSub=mean(Subthreshold_silent(rois{r},t),1)-movmedian(mean(Subthreshold_silent(rois{r},t),1,'omitnan'),smooth_window,'omitnan');
    dSub=movmean(dSub,5);
    [~, upward_trace(r,t)]=detect_transient(dSub,[up_thres 1],zeros(1,length(t)));
    [~, downward_trace(r,t)]=detect_transient(-dSub,[dw_thres 1],zeros(1,length(t)));

    figure(12); clf; tiledlayout(3,3)
    nexttile([1 1])
    histogram(upward_pk{r,1},binEdge); hold all
    histogram(downward_pk{r,1},binEdge);
    title([title_str{r} ', during blue on'])
    legend({'Upward','Downward'})

    nexttile([1 1])
    histogram(upward_pk{r,2},binEdge); hold all
    histogram(downward_pk{r,2},binEdge);
    title([title_str{r} ', during blue off'])
    legend({'Upward','Downward'})

    nexttile([1 1])
    plot(MedianVRposition',squeeze(MedianF(:,:,r))','.')



%%
Lap_excitation=[];  Lap_inhibition=[];
upward_trace(:,Result.Blue>0)=NaN; downward_trace(:,Result.Blue==0)=NaN;
for r=1:3
    Lap_excitation(:,:,r)=PlaceTrigger_average(upward_trace(r,:)>0,nBin,Result.VR,0,115,'rate');
    Lap_inhibition(:,:,r)=PlaceTrigger_average(downward_trace(r,:)>0,nBin,Result.VR,0,115,'rate');
end
figure(13); clf; tiledlayout(3,3)
for r=1:3
    nexttile(r,[1 1])
    imagesc(Lap_excitation(:,:,r))
    title([title_str{r} ,' excitation'])

    nexttile(3+r,[1 1])
    imagesc(Lap_inhibition(:,:,r))
    title([title_str{r} ,' inhibition'])

    nexttile(6+r,[1 1])
    plot([1:nBin]*2/nBin,mean(ringmovMean(Lap_excitation(:,:,r),5),1,'omitnan')); hold all
    plot([1:nBin]*2/nBin,mean(ringmovMean(Lap_inhibition(:,:,r),5),1,'omitnan'));
    yyaxis right;
    plot([1:nBin]*2/nBin,mean(ringmovMean(Lap_FR*1000,3),1,'omitnan'),'k');
    legend({'Excitation','Inhibition','Firing rate'})
    xlabel('VR position (m)')
end
%%
figure(14); clf; tiledlayout(6,7)
NormalizedTrace_filt=movmean(NormalizedTrace_dirt,5,2,'omitnan');
for r=1:3
    subMov=NormalizedTrace_filt(Result.dist_order(noi_dist),find(upward_trace(r,:)>0));
    subMov(isnan(subMov))=0;
    covMat=subMov*subMov';

    [V_ext, D_ext] = eig(covMat);
    D_ext = diag(D_ext);
    D_ext = D_ext(end:-1:1);
    V_ext = V_ext(:,end:-1:1);
    vSign = sign(max(V_ext) - max(-V_ext));  % make the largest value always positive
    V_ext = V_ext.*vSign;

    subMov=NormalizedTrace_filt(Result.dist_order(noi_dist),find(downward_trace(r,:)>0));
    subMov(isnan(subMov))=0;
    covMat=subMov*subMov';
    [V_inh, D_inh] = eig(covMat);
    D_inh = diag(D_inh);
    D_inh = D_inh(end:-1:1);
    V_inh = V_inh(:,end:-1:1);
    vSign = sign(max(V_inh) - max(-V_inh));  % make the largest value always positive
    V_inh = V_inh.*vSign;

    for n=1:7
        nexttile([1 1])
        imshow2(max(double(Result.ftprnt(:,:,Result.dist_order(noi_dist))>0).*reshape(rescale(V_ext(:,n)),1,1,[]),[],3),[])
        colormap('turbo')
        title([title_str{r} ', ext., PC' num2str(n) ', ' num2str(D_ext(n)/sum(D_ext)*100,2) '%'])
    end

    for n=1:7
        nexttile([1 1])
        imshow2(max(double(Result.ftprnt(:,:,Result.dist_order(noi_dist))>0).*reshape(rescale(V_inh(:,n)),1,1,[]),[],3),[])
        colormap('turbo')
        title([title_str{r} ', inh., PC' num2str(n) ', ' num2str(D_inh(n)/sum(D_inh)*100,2) '%'])
    end
end

%% Find modes
noi=setdiff([1:nROI],[BadROI{f}]); noi_dist=ismember(Result.dist_order,noi);
nTauClass={[-40:-2],[-40:-2],[-40:-2]};
nTau_silent=[-1:25]; time_sub=[]; frac_th=0.95;
Spclass_trace=zeros(1,size(Result.traces,2));
ss_time=find(Result.SpClass(1,:));
ISI_ss=ss_time(2:end)-ss_time(1:end-1); bwBrst=bwlabel(ISI_ss<15);
[~,ind]=unique(bwBrst); ind=ind(2:end);
brst_time=ss_time(ind);
ss_time(find(bwBrst>0)+1)=[];
Spclass_trace(ss_time)=1; %SS
Spclass_trace(brst_time)=2; %Brst
Spclass_trace(find(Result.SpClass(2,:)==1))=3; %CS

for c=1:3
    time_sub{c}=tovec(find(Spclass_trace==c)'+nTauClass{c});
end

Silent_trace=NormalizedTrace_dirt(Result.dist_order(noi_dist),:);
Silent_trace(:,unique([tovec(find(Result.SpClass(1,:))'+nTau_silent)' find(Result.CStrace>0)]))=NaN; %Remove peri-spike times
Silent_traceOn=Silent_trace(:,Result.Blue>0); Silent_traceOn(:,sum(isnan(Silent_traceOn),1)>0)=[];
Silent_traceOff=Silent_trace(:,Result.Blue==0); Silent_traceOff(:,sum(isnan(Silent_traceOff),1)>0)=[];
[V_silentOn D_silentOn]=get_eigvector(Silent_traceOn); n_silentOn=sum(cumsum(D_silentOn)/sum(D_silentOn)>frac_th==0);
[V_silentOff D_silentOff]=get_eigvector(Silent_traceOff); n_silentOff=sum(cumsum(D_silentOff)/sum(D_silentOff)>frac_th==0);

prespike_trace=NormalizedTrace_dirt(Result.dist_order(noi_dist),:);
prespike_trace(:,setdiff([1:size(Result.traces,2)],unique(cell2mat(time_sub'))))=NaN; %Take pre-spike dynamics
prespike_traceOn=prespike_trace(:,Result.Blue>0); prespike_traceOn(:,sum(isnan(prespike_traceOn),1)>0)=[];
prespike_traceOff=prespike_trace(:,Result.Blue==0); prespike_traceOff(:,sum(isnan(prespike_traceOff),1)>0)=[];
[V_preSpOn D_preSpOn]=get_eigvector(prespike_traceOn); n_preSpOn=sum(cumsum(D_preSpOn)/sum(D_preSpOn)>frac_th==0);
[V_preSpOff D_preSpOff]=get_eigvector(prespike_traceOff); n_preSpOff=sum(cumsum(D_preSpOff)/sum(D_preSpOff)>frac_th==0);

figure(14); clf; tiledlayout(4,6)
for n=1:6 %Silent during blue on
    nexttile([1 1])
    imshow2(max(double(Result.ftprnt(:,:,Result.dist_order(noi_dist))>0).*reshape(rescale(V_silentOn(:,n)),1,1,[]),[],3),[])
    colormap('turbo')
    if n==1;
        ylabel(['During silent' newline 'Blue On']);
    end
    title(['PC' num2str(n) ', ' num2str(D_silentOn(n)/sum(D_silentOn)*100,2) '%'])
end

for n=1:6 %Silent during blue off
    nexttile([1 1])
    imshow2(max(double(Result.ftprnt(:,:,Result.dist_order(noi_dist))>0).*reshape(rescale(V_silentOff(:,n)),1,1,[]),[],3),[])
    colormap('turbo')
    if n==1;
        ylabel(['During silent' newline 'Blue Off']);
    end
    title(['PC' num2str(n) ', ' num2str(D_silentOff(n)/sum(D_silentOff)*100,2) '%'])
end

for n=1:6 %Pre-spike during blue on
    nexttile([1 1])
    imshow2(max(double(Result.ftprnt(:,:,Result.dist_order(noi_dist))>0).*reshape(rescale(V_preSpOn(:,n)),1,1,[]),[],3),[])
    colormap('turbo')
    if n==1;
        ylabel(['Pre-spike' newline 'Blue On']);
    end
    title(['PC' num2str(n) ', ' num2str(D_preSpOn(n)/sum(D_preSpOn)*100,2) '%'])
end

for n=1:6 %Pre-spike during blue off
    nexttile([1 1])
    imshow2(max(double(Result.ftprnt(:,:,Result.dist_order(noi_dist))>0).*reshape(rescale(V_preSpOff(:,n)),1,1,[]),[],3),[])
    colormap('turbo')
    if n==1;
        ylabel(['Pre-spike' newline 'Blue Off']);
    end
    title(['PC' num2str(n) ', ' num2str(D_preSpOff(n)/sum(D_preSpOff)*100,2) '%'])
end

Vs_reduce=V_preSpOn(:,1:7);
% Vs=[V_silentOn(:,1:n_silentOn) V_silentOff(:,1:n_silentOff) V_preSpOn(:,1:n_preSpOn) V_preSpOff(:,1:n_preSpOff)];
% Ds=[D_silentOn(1:n_silentOn)'/sum(D_silentOn) D_silentOff(1:n_silentOff)'/sum(D_silentOff) ...
%     D_preSpOn(1:n_preSpOn)'/sum(D_preSpOn) D_preSpOff(1:n_preSpOff)'/sum(D_preSpOff)];
% Vs=rescale2(Vs.*sign(skewness(Vs,0,1)),1);
% distances = pdist(Vs', 'euclidean');  % Euclidean distance
% Z = linkage(distances, 'average');
% num_clusters = 7;
% cluster_indices = cluster(Z, 'maxclust', num_clusters);
% Vs_reduce = zeros(size(Vs, 1), num_clusters); cluster_weight=zeros(1,num_clusters);
% weight_indices=cluster_indices.*Ds';
% for i = 1:num_clusters
%     %Vs_reduce(:, i) = mean(Vs(:, find(cluster_indices == i)),2);
%     Vs_reduce(:, i) = Vs(:, find(cluster_indices == i,1));
%     cluster_weight(i)=mean(weight_indices(find(cluster_indices == i)));
% end
% 
% [~,sort_cluster]=sort(cluster_weight,'descend');
% Vs_reduce=Vs_reduce(:,sort_cluster);
% % Alternatively, you can plot a dendrogram to visualize the clustering
% %imagesc(rescale2(Vs_reduce(:,sort_cluster),1))
% figure(15); clf; tiledlayout(1,num_clusters);
% for n=1:num_clusters %Silent during blue off
%     nexttile([1 1])
%     imshow2(max(double(Result.ftprnt(:,:,Result.dist_order(noi_dist))>0).*reshape(rescale(Vs_reduce(:,n)),1,1,[]),[],3),[])
%     colormap('turbo')
%     title(['Cluster' num2str(n)])
% end
%% Mode trace, statistics
noi_rois=zeros(1,nROI); for r=1:3; noi_rois(rois{r})=r; end
noi_rois=noi_rois(Result.dist_order(noi_dist));
Vs_simple=zeros(length(noi),3); 
Vs_simple([find(noi_rois==1) find(noi_rois==2)],1)=1; % first mode, basal and apical -> [1 1]
Vs_simple([find(noi_rois==1)],2)=1; Vs_simple([find(noi_rois==2)],2)=-1; % second mode, basal and apical -> [1 -1]
Vs_simple([find(noi_rois==3)],3)=1; % third mode, peri-soma

%ModeTrace=(NormalizedTrace_dirt(Result.dist_order(noi_dist),:)'*(Vs_reduce))';
ModeTrace=(NormalizedTrace_dirt(Result.dist_order(noi_dist),:)'*(Vs_simple))';
nMode=size(ModeTrace,1);

STA_CSModemat=reshape(ModeTrace(1:nMode,CS_s'+nTau{2}),nMode,[],length(nTau{2}));
STA_SSModemat=reshape(ModeTrace(1:nMode,bAP_s'+nTau{1}),nMode,[],length(nTau{1}));
STA_dSPModemat=reshape(ModeTrace(1:nMode,dSP_s'+nTau{3}),nMode,[],length(nTau{3}));
ModeTrace_prespike=NaN(size(ModeTrace));
ModeTrace_prespike(:,bAP_s'+nTau{1}(1:-nTau{1}(1)-1))=ModeTrace(:,bAP_s'+nTau{1}(1:-nTau{1}(1)-1));
ModeTrace_prespike(:,CS_s'+nTau{2}(1:-nTau{2}(1)-1))=ModeTrace(:,CS_s'+nTau{2}(1:-nTau{2}(1)-1));

ModeTrace_slient=ModeTrace;
t_sp=find(Result.SpClass(1,:));
sp_na=sum((t_sp'+nTau_silent)<0 | (t_sp'+nTau_silent)>size(Result.traces,2),2)==0; t_sp(~sp_na)=[];
ModeTrace_slient(:,unique([tovec(t_sp'+nTau_silent)' find(Result.CStrace>0)]))=NaN; %Remove peri-spike times

runtime=Result.VR(end,:)>0.005 & Result.Blue ==0; resttime=Result.VR(end,:)<0.001 & Result.Blue ==0;

    STA_SSModemat_run=reshape(ModeTrace(1:nMode,bAP_s(find(runtime(bAP_s)))'+nTau{1}),nMode,[],length(nTau{1}));
    STA_SSModemat_rest=reshape(ModeTrace(1:nMode,bAP_s(find(resttime(bAP_s)))'+nTau{1}),nMode,[],length(nTau{1}));
    STA_CSModemat_run=reshape(ModeTrace(1:nMode,CS_s(find(runtime(CS_s)))'+nTau{2}),nMode,[],length(nTau{2}));
    STA_CSModemat_rest=reshape(ModeTrace(1:nMode,CS_s(find(resttime(CS_s)))'+nTau{2}),nMode,[],length(nTau{2}));
    STA_dSPModemat_run=reshape(ModeTrace(1:nMode,dSP_s(find(runtime(dSP_s)))'+nTau{3}),nMode,[],length(nTau{3}));
    STA_dSPModemat_rest=reshape(ModeTrace(1:nMode,dSP_s(find(resttime(dSP_s)))'+nTau{3}),nMode,[],length(nTau{3}));

    STA_SSModemat_on=reshape(ModeTrace(1:nMode,bAP_s(find(Result.Blue(bAP_s)>0))'+nTau{1}),nMode,[],length(nTau{1}));
    STA_SSModemat_off=reshape(ModeTrace(1:nMode,bAP_s(find(Result.Blue(bAP_s)==0))'+nTau{1}),nMode,[],length(nTau{1}));
    STA_CSModemat_on=reshape(ModeTrace(1:nMode,CS_s(find(Result.Blue(CS_s)>0))'+nTau{2}),nMode,[],length(nTau{2}));
    STA_CSModemat_off=reshape(ModeTrace(1:nMode,CS_s(find(Result.Blue(CS_s)==0))'+nTau{2}),nMode,[],length(nTau{2}));
    STA_dSPModemat_on=reshape(ModeTrace(1:nMode,dSP_s(find(Result.Blue(dSP_s)>0))'+nTau{3}),nMode,[],length(nTau{3}));
    STA_dSPModemat_off=reshape(ModeTrace(1:nMode,dSP_s(find(Result.Blue(dSP_s)==0))'+nTau{3}),nMode,[],length(nTau{3}));

    % h=histogram(ModeTrace(n,runtime),100,'Normalization','probability'); hold all
    % histogram(ModeTrace(n,resttime),h.BinEdges,'Normalization','probability'); hold all

    Lap_Mode_prespike = PlaceTrigger_average(movmean(ModeTrace_prespike,7,2,'omitnan'),150,Result.VR,0,115); %total trace
    Lap_Mode_silent = PlaceTrigger_average(movmean(ModeTrace_slient,20,2,'omitnan'),150,Result.VR,0,115); %total trace
    

for n=1:nMode

    figure(25+n); clf; cmap=distinguishable_colors(6); ax1=[];
cax=[-2 2];%*0.01;
tiledlayout(3,4);

ax_mode=nexttile([1 1]);
show_mode=max(double(Result.ftprnt(:,:,Result.dist_order(noi_dist))>0).*reshape(rescale(Vs_simple(:,n)),1,1,[]),[],3);
show_mode(max(Result.ftprnt(:,:,Result.dist_order(noi_dist)),[],3)==0)=0.5;
imshow2(show_mode,[])
colormap(ax_mode,gen_colormap([0 0.5 1; 1 1 1; 1 0 0]))

    ax1=[];
    ax1=[ax1 nexttile([1 1])];
plot(nTau{1},squeeze(mean(STA_SSModemat(n,:,:),2,'omitnan')),'color',[0 0 0]); hold all
plot(nTau{2},squeeze(mean(STA_CSModemat(n,:,:),2,'omitnan')),'color',[1 0 0.2]); hold all
plot(nTau{3},squeeze(mean(STA_dSPModemat(n,:,:),2,'omitnan')),'color',[0 0.6 1]); hold all
legend({['SS, N=' num2str(size(STA_SSModemat,2))],['CS, N=' num2str(size(STA_CSModemat,2))],['dSP, N=' num2str(size(STA_dSPModemat,2))]})
title([num2str(n) 'th Mode'])

    ax1=[ax1 nexttile([1 1])];
plot(nTau{1},squeeze(mean(STA_SSModemat_run(n,:,:),2,'omitnan')),'color',[0 0 0]); hold all
plot(nTau{2},squeeze(mean(STA_CSModemat_run(n,:,:),2,'omitnan')),'color',[1 0 0.2]); hold all
plot(nTau{3},squeeze(mean(STA_dSPModemat_run(n,:,:),2,'omitnan')),'color',[0 0.6 1]); hold all
legend({['SS, N=' num2str(size(STA_SSModemat_run,2))],['CS, N=' num2str(size(STA_CSModemat_run,2))],['dSP, N=' num2str(size(STA_dSPModemat_run,2))]})
title([num2str(n) 'th, Running'])

 ax1=[ax1 nexttile([1 1])];
plot(nTau{1},squeeze(mean(STA_SSModemat_rest(n,:,:),2,'omitnan')),'color',[0 0 0]); hold all
plot(nTau{2},squeeze(mean(STA_CSModemat_rest(n,:,:),2,'omitnan')),'color',[1 0 0.2]); hold all
plot(nTau{3},squeeze(mean(STA_dSPModemat_rest(n,:,:),2,'omitnan')),'color',[0 0.6 1]); hold all
legend({['SS, N=' num2str(size(STA_SSModemat_rest,2))],['CS, N=' num2str(size(STA_CSModemat_rest,2))],['dSP, N=' num2str(size(STA_dSPModemat_rest,2))]})
title([num2str(n) 'th, Resting'])

 ax1=[ax1 nexttile([1 1])];
plot(nTau{1},squeeze(mean(STA_SSModemat_on(n,:,:),2,'omitnan')),'color',[0 0 0]); hold all
plot(nTau{2},squeeze(mean(STA_CSModemat_on(n,:,:),2,'omitnan')),'color',[1 0 0.2]); hold all
plot(nTau{3},squeeze(mean(STA_dSPModemat_on(n,:,:),2,'omitnan')),'color',[0 0.6 1]); hold all
legend({['SS, N=' num2str(size(STA_SSModemat_on,2))],['CS, N=' num2str(size(STA_CSModemat_on,2))],['dSP, N=' num2str(size(STA_dSPModemat_on,2))]})
title([num2str(n) 'th, Blue On'])

 ax1=[ax1 nexttile([1 1])];
plot(nTau{1},squeeze(mean(STA_SSModemat_off(n,:,:),2,'omitnan')),'color',[0 0 0]); hold all
plot(nTau{2},squeeze(mean(STA_CSModemat_off(n,:,:),2,'omitnan')),'color',[1 0 0.2]); hold all
plot(nTau{3},squeeze(mean(STA_dSPModemat_off(n,:,:),2,'omitnan')),'color',[0 0.6 1]); hold all
legend({['SS, N=' num2str(size(STA_SSModemat_off,2))],['CS, N=' num2str(size(STA_CSModemat_off,2))],['dSP, N=' num2str(size(STA_dSPModemat_off,2))]})
title([num2str(n) 'th, Blue Off'])

nexttile([1 1]);
    h=histogram(ModeTrace(n,Result.Blue==0),100,'Normalization','probability'); hold all
    histogram(ModeTrace(n,Result.Blue>0),h.BinEdges,'Normalization','probability'); hold all
    legend({'Blue Off','Blue On'})
linkaxes([ax1],'xy')

ax_pf(1)=nexttile([1 2]);
Lap_Mode_silent_mean=mean(Lap_Mode_silent(:,:,n),1,'omitnan');
show_im=Lap_Mode_silent(:,:,n);
show_im(isnan(show_im))=median(show_im(:),'omitnan');
imagesc(show_im)
title(['Mode #' num2str(n) ',Silent']); colorbar;
xlabel('Position (bin)')
ylabel('Lap')

ax_pf(2)=nexttile([1 2]);
Lap_Mode_prespike_mean=mean(Lap_Mode_prespike(:,:,n),1,'omitnan');
show_im=Lap_Mode_prespike(:,:,n);
show_im(isnan(show_im))=median(show_im(:),'omitnan');
imagesc(show_im)
arrayfun(@(x) colormap(x,gen_colormap([0 0.5 1; 1 1 1; 1 0 0])),ax_pf)
title(['Mode #' num2str(n) ', pre-spike']); colorbar;
xlabel('Position (bin)')
ylabel('Lap')

nexttile([1 1]);
plot([1:nBin]*2/nBin,movmean(Lap_Mode_silent_mean,7,'omitnan')); hold all 
yyaxis right
plot([1:nBin]*2/nBin,movmean(Lap_Mode_prespike_mean,7,'omitnan'),'LineStyle','-'); hold all 
xlabel('Position (m)')
legend({'Silent','Pre-spike'})

end

%%
figure(11); clf; tiledlayout(1,3);
nexttile([1 1]);
histogram(mean(Subthreshold_silent(rois{3},Result.Blue==0),1),h.BinEdges,'Normalization','probability'); hold all
histogram(mean(Subthreshold_silent(rois{3},Result.Blue>0),1),h.BinEdges,'Normalization','probability'); hold all
legend({'Blue off','Blue on'});
title('Peri-Soma subthreshold')

nexttile([1 1]);
h=histogram(mean(Subthreshold_silent(rois{1},Result.Blue==0),1),100,'Normalization','probability'); hold all
histogram(mean(Subthreshold_silent(rois{1},Result.Blue>0),1),h.BinEdges,'Normalization','probability'); hold all
legend({'Blue off','Blue on'});
title('Basal subthreshold')

nexttile([1 1]);
histogram(mean(Subthreshold_silent(rois{2},Result.Blue==0),1),h.BinEdges,'Normalization','probability'); hold all
histogram(mean(Subthreshold_silent(rois{2},Result.Blue>0),1),h.BinEdges,'Normalization','probability'); hold all
legend({'Blue off','Blue on'});
title('Apical subthreshold')






figure(15); clf;
BinEdge=[-10:0.1:10];
tiledlayout(6,2)
for b=1:6
    nexttile;
    h=histogram(mean(Subthreshold_silent(rois{1},find(bwBlue==b)+2000),1),BinEdge,'Normalization','probability'); hold all
    histogram(mean(Subthreshold_silent(rois{1},bwBlue==b),1),BinEdge,'Normalization','probability'); hold all

    nexttile;
    h=histogram(mean(Subthreshold_silent(rois{2},find(bwBlue==b)+2000),1),BinEdge,'Normalization','probability'); hold all
    histogram(mean(Subthreshold_silent(rois{2},bwBlue==b),1),BinEdge,'Normalization','probability'); hold all
end

%% dSpike footprint
f=21; bound=6;
nROI=size(Result.normTraces,1);
nTau={[-30:20],[-70:100],[-30:20]}; %SS, CS, dSP
spclass_ind=1;
%load aligned movie somatic spike
alignmovlist=dir(fullfile(fpath{f},[alignedMovFN{spclass_ind} '*.tiff']));
AlignMov=[];
for l=1:length(alignmovlist)
    l
    AlignMov=cat(3,AlignMov,readtiff(fullfile(fpath{f},alignmovlist(l).name)));
end
sz_align=size(AlignMov);
AlignMov=double(reshape(AlignMov,sz_align(1),sz_align(2),length(nTau{spclass_ind}),[]));
AlignMov=AlignMov-mean(AlignMov(:,:,1:10,:),3);

valid_sp=find(Result.StackedSpike{spclass_ind}(1,:)>0);
dSPikeMat=max(reshape(Result.spike(:,Result.StackedSpike{spclass_ind}(2,valid_sp)'+[-1:1]),nROI,[],3),[],3);
dSpike_ROI=find(sum(dSPikeMat,2)>0);
F0=imgaussfilt(Result.ref_im(bound:end-bound,bound:end-bound),3);
AlignMov_dFF=AlignMov(:,:,:,valid_sp)./F0;
g=1; dSpikeROImov=cell(1,nROI);
for r=dSpike_ROI'
    dSpikeROImov{r}=AlignMov_dFF(:,:,:,find(dSPikeMat(r,:)));
end

noi_dist=[15 17];
catMov=cell2mat(cellfun(@(x) mat2gray(mean(x,4)),dSpikeROImov(noi_dist),'UniformOutput',false)');
moviefixsc(imgaussfilt(catMov,2))

STA_dSPmat_ROI=[]; ax1=[];
dSp_time=Result.StackedSpike{spclass_ind}(2,valid_sp);
figure; tiledlayout(1,2)
for r=noi_dist%dSpike_ROI'
    STA_dSPmat_ROI(:,:,r)=squeeze(mean(reshape(NormalizedTrace(Result.dist_order,dSp_time(find(dSPikeMat(r,:)))'+nTau{3}),nROI,[],length(nTau{3})),2));
    ax1=[ax1 nexttile([1 1])];
    imagesc(STA_dSPmat_ROI(:,:,r))
    title(['ROI #' num2str(r) ', N = ' num2str(length(find(dSPikeMat(r,:))))])
end
linkaxes(ax1,'xy')

figure; show_footprnt(Result.ftprnt(:,:,:),Result.ref_im)