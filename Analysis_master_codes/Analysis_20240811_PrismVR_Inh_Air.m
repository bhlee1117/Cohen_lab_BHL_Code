clear
clc;
[~, ~, raw] = xlsread(['/Users/bhlee1117/Documents/BHL/Matlab_project/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:P27');
cd '/Users/bhlee1117/Documents/BHL/Matlab_project/20240823_Inhbit'
% [~, ~, NeuronsToUse]=xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
%     'PlaceCellData_Arrangement.xlsx'], 'Sheet1', 'L8:M46');
%
% NeuronsToUse=cellfun(@(x) (str2num(num2str(x))),NeuronsToUse,'UniformOutput',false);
ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
basal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
apical_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false);
fpath=raw(:,1)';
StructureData=raw(:,10);
EndFrame=cell2mat(raw(:,13));
ifmotionReject=cell2mat(raw(:,14));
save_figto='/Users/bhlee1117/Documents/BHL/Matlab_project/20240823_Inhbit';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;

%%
f=21;
load('/Users/bhlee1117/Documents/BHL/Matlab_project/Data/Neuron21_data/PC_Result.mat','Result')
%rois={basal_ROI{f},apical_ROI{f}};
nROI=size(Result.normTraces,1);
nTau_bAP=[-20:20];
nTau={[-70:15],[-70:100],[-50:20]}; %SS, CS, dSP
% Isolated Somatic spike
som_spike=find(Result.spike(1,:));
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
for s=find(Result.SpClass(3,:))
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


STA_SS=squeeze(mean(reshape(Result.normTraces(:,bAP_s'+nTau_bAP),nROI,[],length(nTau_bAP)),2));
STA_SS= STA_SS;% - mean(STA_SS(:,1:10),2);
F_ref=mean(STA_SS(:,-nTau_bAP(1)+[7:14]),2);

NormalizedTrace=(Result.normTraces)./F_ref;
%NormalizedTrace_ch=cellfun(@(x) x./F_ref,Result.norm_trace_check,'UniformOutput',false);

STA_CSmat=reshape(NormalizedTrace(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));
STA_SSmat=reshape(NormalizedTrace(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1}));
STA_dSPmat=reshape(NormalizedTrace(:,dSP_s'+nTau{3}),nROI,[],length(nTau{3}));

%STA_CSmat_ch=cellfun(@(x) reshape(x(:,CS_s'+nTau{2}),nROI,[],length(nTau{2})),Result.norm_trace_check,'UniformOutput',false);
%STA_SSmat_ch=cellfun(@(x) reshape(NormalizedTrace(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1})),Result.norm_trace_check,'UniformOutput',false);
%STA_dSPmat_ch=cellfun(@(x) reshape(NormalizedTrace(:,dSP_s'+nTau{3}),nROI,[],length(nTau{3})),Result.norm_trace_check,'UniformOutput',false);

%%
figure(5); clf; cmap=distinguishable_colors(6); ax1=[];
%rois={basal_ROI{f},apical_ROI{f}}; 
cax=[-0.5 2];
tiledlayout(1,3);
ax3=nexttile([1 1]);
imagesc(nTau{1},[1:nROI],squeeze(mean(STA_SSmat(dist_order,:,:),2)),cax)
title('Simple spike')
ax4=nexttile([1 1]);
imagesc(nTau{2},[1:nROI],squeeze(mean(STA_CSmat(dist_order,:,:),2)),cax); %shading interp;
title('Complex spike')
ax5=nexttile([1 1]);
imagesc(nTau{3},[1:nROI],squeeze(mean(STA_dSPmat(dist_order,:,:),2)),cax)
title('Dendritic spike')
colormap(turbo); linkaxes([ax3 ax4 ax5],'xy')

%%

f=20; bound=6;
nROI=size(Result.normTraces,1);
nTau={[-30:20],[-50:100],[-30:20]}; %SS, CS, dSP
spclass_ind=3;
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

noi=[15 17];
catMov=cell2mat(cellfun(@(x) mat2gray(mean(x,4)),dSpikeROImov(noi),'UniformOutput',false)');
moviefixsc(imgaussfilt(catMov,2))
%%
STA_dSPmat_ROI=[]; ax1=[];
dSp_time=Result.StackedSpike{spclass_ind}(2,valid_sp); 
figure; tiledlayout(1,2)
for r=noi%dSpike_ROI'
STA_dSPmat_ROI(:,:,r)=squeeze(mean(reshape(NormalizedTrace(Result.dist_order,dSp_time(find(dSPikeMat(r,:)))'+nTau{3}),nROI,[],length(nTau{3})),2));
ax1=[ax1 nexttile([1 1])];
imagesc(STA_dSPmat_ROI(:,:,r))
title(['ROI #' num2str(r) ', N = ' num2str(length(find(dSPikeMat(r,:))))])
end
linkaxes(ax1,'xy')

figure; show_footprnt(Result.ftprnt(:,:,:),Result.ref_im)
%% Temporal correlation
rois={basal_ROI{f},apical_ROI{f}}; 
nTau_silent=[-50:200];
exclude_frq=[55.5 56.2];
all_spike=find(max(Result.SpClass([1 2 3],:),[],1));

Fs=1000; NormalizedTrace_filt=[];
freq_lowhigh=exclude_frq/(Fs/2);
[b, a] = butter(4, freq_lowhigh, 'stop');
for n=1:nROI
NormalizedTrace_filt(n,:) = filtfilt(b, a, NormalizedTrace(n,:));
end

SilentTrace=NormalizedTrace_filt;
SilentTrace(:,all_spike'+nTau_silent)=0;
SilentTrace(:,Result.motionReject)=0;
SilentTrace=movmean(SilentTrace,5,2,'omitnan');

[SilentTrace_BA]=[mean(SilentTrace(rois{1},:),1); mean(SilentTrace(rois{2},:),1)];
[autocorrBA(:,1)  lagTau]=xcorr(SilentTrace_BA(1,:),5000);
[autocorrBA(:,2)  lagTau]=xcorr(SilentTrace_BA(2,:),5000);
[crosscorrBA  lagTau]=xcorr(SilentTrace_BA(1,:),SilentTrace_BA(2,:),5000);


figure(50); clf;
plot(lagTau,crosscorrBA,'k'); hold all
plot(lagTau,autocorrBA); hold all

%% 
Subthreshold=get_subthreshold(NormalizedTrace,all_spike,5,15);
Subthreshold_silent=Subthreshold;
Subthreshold_silent(:,all_spike'+[-5:10])=NaN;
Subthreshold_silent(:,Result.motionReject)=NaN;
bwBlue=bwlabel(Result.Blue>0);
figure(10); clf;
nexttile;
h=histogram(mean(Subthreshold_silent(rois{1},Result.Blue==0),1),100,'Normalization','probability'); hold all
histogram(mean(Subthreshold_silent(rois{1},Result.Blue>0),1),h.BinEdges,'Normalization','probability'); hold all

nexttile;
histogram(mean(Subthreshold_silent(rois{2},Result.Blue==0),1),h.BinEdges,'Normalization','probability'); hold all
histogram(mean(Subthreshold_silent(rois{2},Result.Blue>0),1),h.BinEdges,'Normalization','probability'); hold all

figure(15); clf;
BinEdge=[-10:0.1:10];
tiledlayout(5,2)
for b=1:5
nexttile;
h=histogram(mean(Subthreshold_silent(rois{1},find(bwBlue==b)+2000),1),BinEdge,'Normalization','probability'); hold all
histogram(mean(Subthreshold_silent(rois{1},bwBlue==b),1),BinEdge,'Normalization','probability'); hold all    
    
nexttile;
h=histogram(mean(Subthreshold_silent(rois{2},find(bwBlue==b)+2000),1),BinEdge,'Normalization','probability'); hold all
histogram(mean(Subthreshold_silent(rois{2},bwBlue==b),1),BinEdge,'Normalization','probability'); hold all    
end