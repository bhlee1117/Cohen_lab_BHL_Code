
clear
clc;
cd '/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:N164');

save_to='/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);
BadROI=cellfun(@(x) (str2num(num2str(x))),raw(:,13),'UniformOutput',false);

%%
i=[72];%length(fpath)
load([fpath{i} '/Result.mat'])

coord_1d=dim_reduce(get_coord(Result.ftprnt));
[~, Result.dist_order]=sort(coord_1d,'descend');
noi=setdiff([1:size(Result.traces,1)],BadROI{i});
noi_dist=ismember(Result.dist_order,noi);

nROI=size(Result.traces,1);

Normtrace=-Result.traces./prctile(-Result.traces,30,2);
%Normtrace=Normtrace-movprc(Normtrace,2000,20,2);
Normtrace=(Normtrace-prctile(Normtrace,5,2))./(prctile(Normtrace,99,2)-prctile(Normtrace,5,2));

% filt_frq=[150 200]; Normtrace_filt=[];
% freq_lowhigh=filt_frq/(1000/2);
% [b4, a4] = butter(4, freq_lowhigh, 'stop');
% for n=1:size(Normtrace,1)
% Normtrace_filt(n,:) = filtfilt(b4, a4, Normtrace(n,:));
% end

figure(101); clf;
tiledlayout(5,1)
%imagesc(Normtrace)
% nexttile([1 1])
% imagesc(Result.ref_im)
ax1=nexttile([2 1]);
imagesc(Normtrace(Result.dist_order(noi_dist),:),[-0.1 1])
colormap(turbo); axis off

ax1=[ax1 nexttile([2 1])];
n=[2 3 4]; plot(rescale(mean(Normtrace(n,:),1)),'k')
hold all
n=[5 7]; plot(rescale(mean(Normtrace(n,:),1)),'r')
legend({'Proximal','Distal'}); axis off

ax2=nexttile([1 1]);
plot(Result.Blue); axis off

linkaxes([ax1 ax2],'x')
xlim([2900 3600])

%%

SomRP_path=[15 22 54 72 86 92 4];
DenRP_path=[13 23 55 74 88 91 3];
path_cat=[SomRP_path; DenRP_path];
rheo_bin=[0:1.3:10];

figure(28); clf; cmap=distinguishable_colors(6); bin_width=1.3;
tiledlayout(3,2);

filterfreq=[55 70];


g=1; dFF=[];
for f=5:6
load(fullfile(fpath{DenRP_path(f)},'Result'));
Ref_F=tovec(Result.ref_im)'*tovec(Result.ftprnt);
dFF(g)=range(Result.traces(1,:))/Ref_F(1);    
tr_ddStim=rescale(Result.normTraces(1,:));
tr_csddStim=Result.CStrace;
if ~isempty(filterfreq)
    freq_lowhigh=filterfreq/(1000/2);
    [b, a] = butter(4, freq_lowhigh, 'stop');
        tr_ddStim = filtfilt(b, a, tr_ddStim')';
end

ax1=nexttile([1 2]);
plot(tr_ddStim+g,'color',[0 0 0]); hold all
show_cs=tr_ddStim.*tr_csddStim; show_cs(tr_csddStim==0)=NaN;
plot(show_cs+g,'color',[0.85 0 0.25]); hold all
axis off
g=g+1;
end
ax2=nexttile([1 2]);
plot(Result.Blue)
axis off
linkaxes([ax1 ax2],'x')


%% STA, SS, CS, dSP

clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:R31');

ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
basal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
apical_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false);
fpath=raw(:,1)';
StructureData=raw(:,10);
BadROI=cellfun(@(x) (str2num(num2str(x))),raw(:,15),'UniformOutput',false);
EndFrame=cell2mat(raw(:,13));
ifmotionReject=cell2mat(raw(:,14));
ifdirtRemoval=cell2mat(raw(:,16));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
set(0,'DefaultFigureWindowStyle','docked')

f=26; load(fullfile(fpath{f},'PC_Result.mat'),'Result')
cd(fpath{f})
rois={basal_ROI{f},apical_ROI{f},ref_ROI{f}};
nROI=size(Result.normTraces,1);
nTau_bAP=[-15:15];
nTau={[-140:50],[-140:50],[-10:10],[-140:50]}; %SS, CS, dSP, Brst
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
%dSP_list=dSP_list(isBlueoff & ~isBasal);
for s=dSP_list
    isnearby=sum(ismember(s+nTau{3},som_spike))>1;
    isnearbyCS=sum(ismember(s+nTau{3},find(Result.CStrace)))>1;
    isnearbydS=sum(ismember(s+nTau{3},som_spike))>1;
    ispartCS=tr_trace(s)>0;
    %if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS & ~isnearbydS & sum(Result.dirtTrace(:,s+nTau{3}),[1 2])==0
    if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS & ~isnearbydS   
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

STA_SS=squeeze(mean(reshape(Result.normTraces(:,bAP_ref'+nTau_bAP),nROI,[],length(nTau_bAP)),2));
STA_SS= STA_SS - mean(STA_SS(:,1:5),2);
F_ref=mean(STA_SS(:,-nTau_bAP(1)+[4:7]),2);
%F_ref=(tovec(imgaussfilt(Result.ref_im,2))'*tovec(Result.ftprnt)/Result.SpikeHeight_fit(1))';

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
%NormalizedTrace_dirt(Result.dirtTrace>0)=NaN;
%NormalizedTrace_dirt(:,17000:18000)=NaN;
NormalizedTrace_dirt(:,Result.motionReject)=NaN;
NormalizedTrace_ch=cellfun(@(x) x./F_ref,Result.norm_trace_check,'UniformOutput',false);
NormalizedTrace_ch{1}(:,Result.motionReject)=NaN; NormalizedTrace_ch{2}(:,Result.motionReject)=NaN;
%NormalizedTrace_ch{1}(Result.dirtTrace>0)=NaN; NormalizedTrace_ch{2}(Result.dirtTrace>0)=NaN;

STA_CSmat=reshape(NormalizedTrace_dirt(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));
STA_SSmat=reshape(NormalizedTrace_dirt(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1}));
STA_dSPmat=reshape(NormalizedTrace_dirt(:,dSP_s'+nTau{3}),nROI,[],length(nTau{3}));
STA_BSmat=reshape(NormalizedTrace_dirt(:,BS_s'+nTau{4}),nROI,[],length(nTau{4}));

STA_CSmat_ch=cellfun(@(x) reshape(x(:,CS_s'+nTau{2}),nROI,[],length(nTau{2})),NormalizedTrace_ch,'UniformOutput',false);
STA_SSmat_ch=cellfun(@(x) reshape(x(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1})),NormalizedTrace_ch,'UniformOutput',false);
STA_dSPmat_ch=cellfun(@(x) reshape(x(:,dSP_s'+nTau{3}),nROI,[],length(nTau{3})),NormalizedTrace_ch,'UniformOutput',false);
STA_BSmat_ch=cellfun(@(x) reshape(x(:,BS_s'+nTau{4}),nROI,[],length(nTau{4})),NormalizedTrace_ch,'UniformOutput',false);

coord_1d=-dim_reduce(get_coord(Result.ftprnt));
coord_1d=coord_1d-coord_1d(1);
[~, Result.dist_order]=sort(coord_1d,'ascend');
%% SS, CS average, f=26;
figure(4); clf; cmap=distinguishable_colors(6); ax1=[];
noi=setdiff([1:nROI],[BadROI{f}]);
noi_dist=ismember(Result.dist_order,noi);
%rois={basal_ROI{f},apical_ROI{f}};
CoordX=[-184 0 293];
cax=[-0.2 2.5];%*0.01;
tiledlayout(1,6);
ax3=nexttile([1 2]);
imagesc(nTau{1},[1:sum(noi_dist)],squeeze(mean(STA_SSmat(Result.dist_order(noi_dist),:,:)-median(STA_SSmat(Result.dist_order(noi_dist),:,1:-nTau{1}(1)),3,'omitnan'),2,'omitnan')),cax)
xlabel('Peri-spike time (ms)')
title(['Average of Simple spike , N = ' num2str(size(STA_SSmat,2))])
set(gca,'YTick',[1 find(Result.dist_order(noi_dist)==1) sum(noi_dist)],'YTickLabel',num2str(CoordX',3))
ax4=nexttile([1 2]);
imagesc(nTau{2},[1:sum(noi_dist)],squeeze(mean(STA_CSmat(Result.dist_order(noi_dist),:,:)-median(STA_CSmat(Result.dist_order(noi_dist),:,1:-nTau{2}(1)),3,'omitnan'),2,'omitnan')),cax); %shading interp;
xlabel('Peri-spike time (ms)')
set(gca,'YTick',[1 find(Result.dist_order(noi_dist)==1) sum(noi_dist)],'YTickLabel',num2str(CoordX',3))
title(['Average of Complex spike , N = ' num2str(size(STA_CSmat,2))])
colormap(turbo)

cax=[-0.2 1];%*0.01;
ax3=nexttile([1 1]);
imagesc(nTau{1}(-nTau{1}(1)+[-30:5]),[1:sum(noi_dist)],squeeze(mean(STA_SSmat(Result.dist_order(noi_dist),:,(-nTau{1}(1)+[-30:5]))-median(STA_SSmat(Result.dist_order(noi_dist),:,1:-nTau{1}(1)),3,'omitnan'),2,'omitnan')),cax)
xlabel('Peri-spike time (ms)')
set(gca,'YTick',[1 find(Result.dist_order(noi_dist)==1) sum(noi_dist)],'YTickLabel',num2str(CoordX',3))
ax4=nexttile([1 1]);
imagesc(nTau{2}(-nTau{2}(1)+[-30:5]),[1:sum(noi_dist)],squeeze(mean(STA_CSmat(Result.dist_order(noi_dist),:,(-nTau{2}(1)+[-30:5]))-median(STA_CSmat(Result.dist_order(noi_dist),:,1:-nTau{2}(1)),3,'omitnan'),2,'omitnan')),cax); %shading interp;
xlabel('Peri-spike time (ms)')
set(gca,'YTick',[1 find(Result.dist_order(noi_dist)==1) sum(noi_dist)],'YTickLabel',num2str(CoordX',3))
colormap(turbo)


%%
f=20; load(fullfile(fpath{f},'PC_Result.mat'),'Result'); nTau_bAP=[-15:15];
bAP_s=[]; nROI=size(Result.normTraces,1);
noi=setdiff([1:nROI],[BadROI{f}]);
noi_dist=ismember(Result.dist_order,noi);
bAP_ref=[];
coord_1d=-dim_reduce(get_coord(Result.ftprnt));
coord_1d=coord_1d-coord_1d(1);
som_ind=find(Result.dist_order(noi_dist)==1);
som_spike=find(Result.spike(1,:));
tr_ref=Result.normTraces(ref_ROI{f},:);
tr_sub=mean(tr_ref,1)-movprc(mean(tr_ref,1),200,20,2);
tr_sub=get_subthreshold(tr_sub,Result.spike(1,:),5,10);
[trans tr_trace]=detect_transient2(tr_sub,[5 1.5],Result.spike(1,:),15);
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

STA_SS=squeeze(mean(reshape(Result.normTraces(:,bAP_ref'+nTau_bAP),nROI,[],length(nTau_bAP)),2));
STA_SS= STA_SS - mean(STA_SS(:,1:5),2);
F_ref=mean(STA_SS(:,-nTau_bAP(1)+[4:7]),2);
NormalizedTrace=(Result.normTraces)./F_ref;
NormalizedTrace_dirt=NormalizedTrace;
NormalizedTrace_dirt(Result.dirtTrace>0)=NaN;
NormalizedTrace_dirt(:,Result.motionReject)=NaN;
rois={basal_ROI{f},apical_ROI{f},ref_ROI{f}};

noi=setdiff([1:nROI],[BadROI{f}]);
noi_dist=ismember(Result.dist_order,noi);
NormalizedTrace_dirt_filt=NormalizedTrace_dirt(:,1:100000);%-movprc(NormalizedTrace_dirt(:,1:100000),1000,20,2);

figure(21); clf; tiledlayout(6,1); cmap=[0 0.4 1;1 0 0;0 0 0];
nexttile([1 1])
imagesc(NormalizedTrace_dirt_filt(Result.dist_order(noi_dist),55700+[-1400:1400]),[-1 3.5])
axis tight off
nexttile([2 1])
for r=1:3
plot(mean(NormalizedTrace_dirt_filt(rois{r},55700+[-1400:1400]),1,'omitnan'),'color',cmap(r,:)); hold all
end
axis tight off
%plot(Result.VR(5,55700+[-1400:1400])/115-4,'color',[1 0.7 0]);
ylim([-4 8])

nexttile([1 1])
imagesc(NormalizedTrace_dirt_filt(Result.dist_order(noi_dist),78400+[-1400:1400]),[-1 3.5]);
axis tight off

colormap(turbo)
nexttile([2 1])
for r=1:3
plot(mean(NormalizedTrace_dirt_filt(rois{r},78400+[-1400:1400]),1,'omitnan'),'color',cmap(r,:)); hold all
end
%plot(Result.VR(5,78400+[-1400:1400])/115-4,'color',[1 0.7 0]);
axis tight off
ylim([-4 8])

%% Representative CS, SS, dSP, Plateau
f=20; load(fullfile(fpath{f},'PC_Result.mat'),'Result'); nTau_bAP=[-15:15];
bAP_s=[]; nROI=size(Result.normTraces,1);
noi=setdiff([1:nROI],[BadROI{f}]);
noi_dist=ismember(Result.dist_order,noi);
bAP_ref=[];
coord_1d=-dim_reduce(get_coord(Result.ftprnt));
coord_1d=coord_1d-coord_1d(1);
som_ind=find(Result.dist_order(noi_dist)==1);
som_spike=find(Result.spike(1,:));
tr_ref=Result.normTraces(ref_ROI{f},:);
tr_sub=mean(tr_ref,1)-movprc(mean(tr_ref,1),200,20,2);
tr_sub=get_subthreshold(tr_sub,Result.spike(1,:),5,10);
[trans tr_trace]=detect_transient2(tr_sub,[5 1.5],Result.spike(1,:),15);
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

STA_SS=squeeze(mean(reshape(Result.normTraces(:,bAP_ref'+nTau_bAP),nROI,[],length(nTau_bAP)),2));
STA_SS= STA_SS - mean(STA_SS(:,1:5),2);
F_ref=mean(STA_SS(:,-nTau_bAP(1)+[4:7]),2);
NormalizedTrace=(Result.normTraces)./F_ref;
NormalizedTrace_dirt=NormalizedTrace;
NormalizedTrace_dirt(Result.dirtTrace>0)=NaN;
NormalizedTrace_dirt(:,Result.motionReject)=NaN;

CoordX=[-140 0 234];
CSpikeIm=NormalizedTrace(Result.dist_order(noi_dist),55601:55800);
dSpikeIm=NormalizedTrace(Result.dist_order(noi_dist),78460+[-50:49]);

f=25; load(fullfile(fpath{f},'PC_Result.mat'),'Result'); nTau_bAP=[-15:15];
bAP_s=[];
nROI=size(Result.normTraces,1);
noi=setdiff([1:nROI],[BadROI{f}]);
noi_dist=ismember(Result.dist_order,noi);
som_spike=find(Result.spike(1,:));
bAP_ref=[];
coord_1d=-dim_reduce(get_coord(Result.ftprnt));
coord_1d=coord_1d-coord_1d(1);
tr_ref=Result.normTraces(ref_ROI{f},:);
tr_sub=mean(tr_ref,1)-movprc(mean(tr_ref,1),200,20,2);
tr_sub=get_subthreshold(tr_sub,Result.spike(1,:),5,10);
[trans tr_trace]=detect_transient2(tr_sub,[5 1.5],Result.spike(1,:),15);
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

STA_SS=squeeze(mean(reshape(Result.normTraces(:,bAP_ref'+nTau_bAP),nROI,[],length(nTau_bAP)),2));
STA_SS= STA_SS - mean(STA_SS(:,1:5),2);
F_ref=mean(STA_SS(:,-nTau_bAP(1)+[4:7]),2);
NormalizedTrace=(Result.normTraces)./F_ref;
NormalizedTrace_dirt=NormalizedTrace;
NormalizedTrace_dirt(Result.dirtTrace>0)=NaN;
NormalizedTrace_dirt(:,Result.motionReject)=NaN;

Plateau1X=coord_1d(Result.dist_order(noi_dist));
Plateau{1}=NormalizedTrace_dirt(Result.dist_order(noi_dist),43181:43480);

f=18; load(fullfile(fpath{f},'PC_Result.mat'),'Result'); nTau_bAP=[-15:15];
bAP_s=[]; nROI=size(Result.normTraces,1);
noi=setdiff([1:nROI],[BadROI{f}]);
noi_dist=ismember(Result.dist_order,noi);
som_spike=find(Result.spike(1,:));
bAP_ref=[];
coord_1d=dim_reduce(get_coord(Result.ftprnt));
coord_1d=coord_1d-coord_1d(1);
tr_ref=Result.normTraces(ref_ROI{f},:);
tr_sub=mean(tr_ref,1)-movprc(mean(tr_ref,1),200,20,2);
tr_sub=get_subthreshold(tr_sub,Result.spike(1,:),5,10);
[trans tr_trace]=detect_transient2(tr_sub,[5 1.5],Result.spike(1,:),15);
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

STA_SS=squeeze(mean(reshape(Result.normTraces(:,bAP_ref'+nTau_bAP),nROI,[],length(nTau_bAP)),2));
STA_SS= STA_SS - mean(STA_SS(:,1:5),2);
F_ref=mean(STA_SS(:,-nTau_bAP(1)+[4:7]),2);
NormalizedTrace=(Result.normTraces)./F_ref;
NormalizedTrace_dirt=NormalizedTrace;
NormalizedTrace_dirt(:,Result.motionReject)=NaN;

Plateau2X=coord_1d(Result.dist_order(noi_dist));
Plateau{2}=NormalizedTrace_dirt(Result.dist_order(noi_dist),25501:25800);

figure(22); clf; cax=[-1 5]; tiledlayout(1,8);
nexttile([1 1])
imagesc(CSpikeIm,cax)
set(gca,'YTick',[1 som_ind size(CSpikeIm)],'YTickLabel',num2str(CoordX',3))
nexttile([1 1])
imagesc(dSpikeIm,cax)
set(gca,'YTick',[1 som_ind size(CSpikeIm)],'YTickLabel',num2str(CoordX',3))
nexttile([1 3])
imagesc([1:300],Plateau1X*0.94,Plateau{1},cax)
nexttile([1 3])
imagesc([1:300],Plateau2X*0.94,Plateau{2},cax)
colormap('turbo')

%%
k=3; % # of clusters

nTau_silent=[-5:100];
all_spike=find(max(Result.SpClass([1 2 3],:),[],1));
sp_na=find(sum((all_spike'+nTau_silent)<1 ,2)>0 | sum((all_spike'+nTau_silent)>EndFrame(f) ,2)>0);
all_spike(sp_na)=[];
SilentTrace_ch=cellfun(@(x) x(Result.dist_order(noi_dist),:),NormalizedTrace_ch,'UniformOutput',false);
SilentTrace_ch{1}(:,all_spike'+nTau_silent)=NaN; SilentTrace_ch{2}(:,all_spike'+nTau_silent)=NaN;

SilentTrace=NormalizedTrace_dirt;
SilentTrace(:,all_spike'+nTau_silent)=NaN;
nonNan_frame=find(sum(isnan(SilentTrace_ch{1}),1)==0);

CorrMat_Silent=[];
for d1=1:sum(noi_dist)
    for d2=1:sum(noi_dist)
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

% Ct = squeeze(CorrMat_Silent);
% Ct(logical(eye(size(Ct)))) = 1;
% distMatrix = sqrt(2*(1 - Ct));
%Z_ref = linkage(squareform(distMatrix), 'average');
distMatrix=pdist(SilentTrace_ch{1}(:,nonNan_frame),"euclidean");
Z_ref = linkage((distMatrix), 'ward');
leafOrder = optimalleaforder(Z_ref,squareform(distMatrix));
Cluster_ref= switchlabel(cluster(Z_ref, 'maxclust', k));
[~, order_cluster]=sort(Cluster_ref,'ascend');

figure(22); clf; ax_corr=[];
tiledlayout(2,2);
cmap=[0 0.1 1;0 1 0;1 0 0];

ax_fprnt=nexttile([1 1]);
cluster_im=max(double(Result.ftprnt(bound:end-bound,bound:end-bound,Result.dist_order(noi_dist))>0).*reshape(Cluster_ref,1,1,[]),[],3);
cluster_stack=[];
for z=1:k
    cluster_stack(:,:,z)=cluster_im==z;
end
cluster_im=cluster_im.*double(dendrite_bin>0.2);
cluster_im_rgb=grs2rgb(cluster_im,colormap(gen_colormap(cmap)),0.7,3);
cluster_im_rgb=cluster_im_rgb.*mat2gray(Result.ref_imSTA)*2.*double(cluster_im>0.2);
show_footprnt_contour(cluster_stack,cluster_im_rgb,cmap); axis equal tight off

ax_corr=nexttile([1 1]);
imagesc(CorrMat_Silent(order_cluster,order_cluster),[0 0.5]); axis tight equal off; hold all
%title('Correlation during silent, reorder by cluster'); 
colormap(ax_corr,'turbo')
scatter(-0.5*ones(sum(noi_dist),1),[1:sum(noi_dist)],40,cmap(Cluster_ref(order_cluster),:),'filled','MarkerEdgeColor','k','Marker','>','MarkerEdgeAlpha',0)

nexttile([1 1]);
dendrogram(Z_ref,0,'ColorThreshold',600); axis off

nexttile([1 1])
for cl=1:k
ind=triu(CorrMat_Silent,1)>0 & ((Cluster_ref==cl)*(Cluster_ref==cl)');
scatter(tovec(interDendDist(ind))*1.17,tovec(CorrMat_Silent(ind)),20,cmap(cl,:),'filled'); hold all
end
legend({'Basal','Peri-soma','Apical'})
xlabel('Pairwise geodesic distance (\mum)')
ylabel('Correlation coefficient')

%%
