clear
clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:P27');

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
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;

%%
f=20; load(fullfile(fpath{f},'PC_Result.mat'),'Result')
cd(fpath{f})
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
all_spike=find(max(Result.SpClass([1 2 3],:),[],1));
SilentTrace=NormalizedTrace;
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
%% STA movie SS, CS

f=20; bound=6;
nTau={[-30:20],[-50:100],[-30:20]}; %SS, CS, dSP
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
AlignMov=AlignMov-mean(AlignMov(:,:,1:6,:),3);

F0=imgaussfilt(Result.ref_im(bound:end-bound,bound:end-bound),3);
AlignMov_dFF=AlignMov./F0;
 
%%
noi=setdiff([1:nROI],[9]);
SS_all=permute(STA_SSmat(dist_order(noi),:,1:-nTau{1}(1)+5),[1 3 2]);
SS_all=tovec(SS_all); SS_all=SS_all-mean(SS_all,2);

CS_all=permute(STA_CSmat(dist_order(noi),:,-nTau{2}(1)-55:-nTau{2}(1)+55),[1 3 2]);
CS_all=tovec(CS_all); CS_all=CS_all-mean(CS_all,2);

covMat = SS_all*SS_all';
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V_ss = V.*vSign;
figure(4); clf;
for i=1:20
    nexttile([1 1])
    imagesc(reshape(V_ss(:,i),length(noi),[]))
    title([num2str(i) ', Fraction: ' num2str(D(i)/sum(D),2)])
end
colormap(turbo)

covMat = CS_all*CS_all';
[V, D] = eig(covMat);
D = diag(D);
D = D(end:-1:1);
V = V(:,end:-1:1);
vSign = sign(max(V) - max(-V));  % make the largest value always positive
V_cs = V.*vSign;
figure(5); clf;
for i=1:20
    nexttile([1 1])
    imagesc(reshape(V_cs(:,i),length(noi),[]))
    title([num2str(i) ', Fraction: ' num2str(D(i)/sum(D),2)])
end
colormap(turbo)

%%
nTau_silent=[-50:200]; Nbin=70;
nTau2={[-300:100],[-300:100],[-1000:100]};
%all_spike=find(max(Result.spike,[],1));
all_spike=find(max(Result.SpClass([1 2],:),[],1));
SilentTrace=NormalizedTrace_ch;
for ch=1:length(SilentTrace)
SilentTrace{ch}(:,all_spike'+nTau_silent)=0;
SilentTrace{ch}(:,Result.motionReject)=0;
SilentTrace{ch}=movmean(SilentTrace{ch},5,2,'omitnan');
%SilentTrace{ch}=SilentTrace{ch}(:,sum(isnan(SilentTrace{ch}))==0);
end

all_spike2=find(max(Result.SpClass([3],:),[],1));
SilentTrace2=NormalizedTrace_ch;
for ch=1:length(SilentTrace)
%SilentTrace2{ch}(:,all_spike2'+nTau_silent)=NaN;
SilentTrace2{ch}(:,Result.motionReject)=0;
SilentTrace2{ch}=movmean(SilentTrace2{ch},5,2,'omitnan');
%SilentTrace2{ch}=SilentTrace2{ch}(:,sum(isnan(SilentTrace2{ch}))==0);
end

silent_basal=mean(SilentTrace{1}(rois{1},:));
silent_apical=mean(SilentTrace{2}(rois{2},:));
all_basal=mean(SilentTrace2{1}(rois{1},:));
all_apical=mean(SilentTrace2{2}(rois{2},:));

figure(7); clf; ax1=[];
Auto_basal_SS=[]; Auto_basal_CS=[]; 
Auto_basal_all=[]; Auto_basal_all=[]; 
Auto_apical_SS=[]; Auto_apical_CS=[];
Cross_bsap_SS=[]; Cross_bsap_CS=[];
tiledlayout(1,3)
ax1=[ax1 nexttile([1 1])];
[Auto_basal t_lag]=xcorr(silent_basal,silent_basal,10000);
[Auto_basal_all t_lag]=xcorr(all_basal,all_basal,10000);
for s=1:size(SS_presub_basal,2)
[Auto_basal_SS(:,s) t_lag_ss]=xcorr(SS_presub_basal(:,s),SS_presub_basal(:,s),'unbiased');
end
for s=1:size(CS_presub_basal,2)
[Auto_basal_CS(:,s) t_lag_cs]=xcorr(CS_presub_basal(:,s),CS_presub_basal(:,s),'unbiased');
end
plot(t_lag,Auto_basal); hold all
plot(t_lag,Auto_basal_all); hold all
% plot(t_lag_ss,mean(Auto_basal_SS,2))
% plot(t_lag_cs,mean(Auto_basal_CS,2))
xlabel('Lag time (ms)')
title('Auto-correlation of subthreshold dynamics of basal dendrite')

ax1=[ax1 nexttile([1 1])];
[Auto_apical t_lag]=xcorr(silent_apical,silent_apical,10000);
[Auto_apical_all t_lag]=xcorr(all_apical,all_apical,10000);
for s=1:size(SS_presub_apical,2)
[Auto_apical_SS(:,s) t_lag_ss]=xcorr(SS_presub_apical(:,s),SS_presub_apical(:,s));
end
for s=1:size(CS_presub_basal,2)
[Auto_apical_CS(:,s) t_lag_cs]=xcorr(CS_presub_apical(:,s),CS_presub_apical(:,s)');
end
plot(t_lag,Auto_apical); hold all
plot(t_lag,Auto_apical_all); hold all
% plot(t_lag_ss,mean(Auto_apical_SS,2))
% plot(t_lag_cs,mean(Auto_apical_CS,2))
xlabel('Lag time (ms)')
title('Auto-correlation of subthreshold dynamics of apical dendrite')

ax1=[ax1 nexttile([1 1])];
[Cross_bsap t_lag]=xcorr(silent_basal,silent_apical,10000);
[Cross_bsap_all t_lag]=xcorr(all_basal,all_apical,10000);
for s=1:size(SS_presub_apical,2)
[Cross_bsap_SS(:,s) t_lag_ss]=xcorr(SS_presub_basal(:,s),SS_presub_apical(:,s));
end
for s=1:size(CS_presub_basal,2)
[Cross_bsap_CS(:,s) t_lag_cs]=xcorr(CS_presub_basal(:,s),CS_presub_apical(:,s));
end
plot(t_lag,Cross_bsap); hold all
plot(t_lag,Cross_bsap_all); hold all
% plot(t_lag_ss,mean(Cross_bsap_SS,2,'omitnan')); hold all
% plot(t_lag_cs,mean(Cross_bsap_CS,2,'omitnan'))
xlabel('Lag time (ms)')
title('Cross-correlation between basal and apical')
legend({'Silent trace','Whole trace'})
%legend({'Pre SS','Pre CS'})
linkaxes(ax1,'x')

%% correlation between dendrites
nTau_silent=[-50:200]; Nbin=70;
all_spike=find(max(Result.spike,[],1));
%all_spike=find(max(Result.SpClass([],:),[],1));
SilentTrace=NormalizedTrace_ch;
for ch=1:length(SilentTrace)
SilentTrace{ch}(:,all_spike'+nTau_silent)=NaN;
SilentTrace{ch}(:,Result.motionReject)=NaN;
SilentTrace{ch}=SilentTrace{ch}(:,sum(isnan(SilentTrace{ch}))==0);
end

Corr_dendrite=[];
for i=1:nROI
    for j=1:nROI
        Corr_dendrite(i,j)=corr(SilentTrace{1}(i,:)',SilentTrace{2}(j,:)');
    end
end
figure(9); clf;
imagesc(Corr_dendrite(dist_order,dist_order))

%%

presp_time=20; bin=1; noi=setdiff([1:nROI],[]);
SS_presub=cellfun(@(x) permute(x(dist_order(noi),:,-nTau{1}(1)-presp_time:-nTau{1}(1)-2),[1 3 2]),STA_SSmat_ch,'UniformOutput',false);
SS_presub=cellfun(@(x) movmean(x,bin,2,'omitnan'),SS_presub,'UniformOutput',false);
for ch=1:2; SS_presub{ch}=SS_presub{ch}(:,[1:bin:end],:); end;

CS_presub=cellfun(@(x) permute(x(dist_order(noi),:,-nTau{2}(1)-presp_time:-nTau{2}(1)-2),[1 3 2]),STA_CSmat_ch,'UniformOutput',false);
CS_presub=cellfun(@(x) movmean(x,bin,2,'omitnan'),CS_presub,'UniformOutput',false);

dSP_presub=cellfun(@(x) permute(x(dist_order(noi),:,-nTau{3}(1)-presp_time:-nTau{3}(1)-2),[1 3 2]),STA_dSPmat_ch,'UniformOutput',false);
dSP_presub=cellfun(@(x) movmean(x,bin,2,'omitnan'),dSP_presub,'UniformOutput',false);

for ch=1:2
    SS_presub{ch}=SS_presub{ch}(:,[1:bin:end],:);
    CS_presub{ch}=CS_presub{ch}(:,[1:bin:end],:); 
    dSP_presub{ch}=dSP_presub{ch}(:,[1:bin:end],:);
end

% CS_presub=permute(STA_CSmat(dist_order(noi),:,-nTau{2}(1)-presp_time:-nTau{2}(1)-2),[1 3 2]);
% CS_presub=movmean(CS_presub,bin,2,'omitnan');
% CS_presub=CS_presub(:,[1:bin:end],:);
% 
% dSP_presub=permute(STA_dSPmat(dist_order(noi),:,-nTau{3}(1)-presp_time:-nTau{3}(1)-2),[1 3 2]);
% dSP_presub=movmean(dSP_presub,bin,2,'omitnan');
% dSP_presub=dSP_presub(:,[1:bin:end],:);

Corr_dendrite_SS=[]; Corr_dendrite_CS=[]; Corr_dendrite_dSP=[];
dist_mat=interDendDist(dist_order,dist_order);
for i=1:nROI
    for j=1:nROI
        for s=1:size(SS_presub{1},3)
        Corr_dendrite_SS(i,j,s)=corr(SS_presub{1}(i,:,s)',SS_presub{2}(j,:,s)');
        end
        for s=1:size(CS_presub{1},3)
        Corr_dendrite_CS(i,j,s)=corr(CS_presub{1}(i,:,s)',CS_presub{2}(j,:,s)');
        end
        for s=1:size(dSP_presub{1},3)
        Corr_dendrite_dSP(i,j,s)=corr(dSP_presub{1}(i,:,s)',dSP_presub{2}(j,:,s)');
        end
    end
end
figure(8); clf; cax=[-0.1 0.7];
cmap=distinguishable_colors(6); 
cmap_dot=cmap+0.4; cmap_dot(cmap_dot>1)=1;
MeancorrMat={mean(Corr_dendrite_SS,3),mean(Corr_dendrite_CS,3),mean(Corr_dendrite_dSP,3)};
tiledlayout(2,3)
nexttile([1 1])
imagesc(mean(Corr_dendrite_SS,3),cax)
title('SS')
nexttile([1 1])
imagesc(mean(Corr_dendrite_CS,3),cax)
title('CS')
nexttile([1 1])
imagesc(mean(Corr_dendrite_dSP,3),cax)
colormap(turbo)
title('dSpike')

nexttile([1 3])
dist_bin=30;
upptri_ind=find(triu(true([length(noi) length(noi)]),1));
x=floor(dist_mat(upptri_ind)/dist_bin);
M=[]; S=[];
for i=unique(x)'
    bin_ind=find(x==i);
    M=[M;[mean(MeancorrMat{1}(upptri_ind(bin_ind))) mean(MeancorrMat{2}(upptri_ind(bin_ind))) mean(MeancorrMat{3}(upptri_ind(bin_ind)))]];
    S=[S;[std(MeancorrMat{1}(upptri_ind(bin_ind))) std(MeancorrMat{2}(upptri_ind(bin_ind))) std(MeancorrMat{3}(upptri_ind(bin_ind)))]];
end
for i=1:3
scatter(dist_mat(upptri_ind),MeancorrMat{i}(upptri_ind),100,cmap_dot(i,:),'.'); hold all
end
hold all
for i=1:3
l(i)=errorbar((unique(x)'+1/2)*dist_bin,M(:,i),S(:,i),'color',cmap(i,:),'LineWidth',2);
end
legend(l(1:3),{'SS','CS','dSP'})
xlabel('Distance (\mum)')
ylabel('Correlation coefficient')
%%
pretime=20;
basalCS=squeeze(mean(STA_CSmat(rois{1},:,-nTau{2}(1)-pretime:-nTau{2}(1)-1)));
apicalCS=squeeze(mean(STA_CSmat(rois{2},:,-nTau{2}(1)-pretime:-nTau{2}(1)-1)));

basalSS=squeeze(mean(STA_SSmat(rois{1},:,-nTau{1}(1)-pretime:-nTau{1}(1)-1)));
apicalSS=squeeze(mean(STA_SSmat(rois{2},:,-nTau{1}(1)-pretime:-nTau{1}(1)-1)));

swingCS=[basalCS(:,end)-mean(mink(basalCS,3,2),2) apicalCS(:,end)-mean(mink(apicalCS,3,2),2)];
swingSS=[basalSS(:,end)-mean(mink(basalSS,3,2),2) apicalSS(:,end)-mean(mink(apicalSS,3,2),2)];

minCS=[mean(mink(basalCS,3,2),2) mean(mink(apicalCS,3,2),2)];
minSS=[mean(mink(basalSS,3,2),2) mean(mink(apicalSS,3,2),2)];

figure(9); clf; ax1=[];
ax1=[ax1 nexttile([1 1])];
scatter(swingSS(:,1),swingSS(:,2),150,'.'); hold all
scatter(swingCS(:,1),swingCS(:,2),150,'.');
xlabel('Basal Swing')
ylabel('Apical Swing')
legend({'SS','CS'});

ax1=[ax1 nexttile([1 1])];
scatter(minSS(:,1),minSS(:,2),150,'.'); hold all
scatter(minCS(:,1),minCS(:,2),150,'.');
xlabel('Min Basal')
ylabel('Min Apical')
legend({'SS','CS'});

ax1=[ax1 nexttile([1 1])];
scatter(basalSS(:,end),apicalSS(:,end),150,'.'); hold all
scatter(basalCS(:,end),apicalCS(:,end),150,'.');
xlabel('Just before spike, Basal')
ylabel('Just before spike, Apical')
legend({'SS','CS'});
linkaxes(ax1,'xy')
%%
CSspikemat=reshape(Result.spike(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));

CSconvmat=conv2(STA_SS',squeeze(CSspikemat(1,8,:)))';
CSconvmat=CSconvmat(dist_order,floor(size(STA_SS,2)/2)+1:end-floor(size(STA_SS,2)/2));
imagesc(squeeze(STA_CSmat(dist_order,8,:))-CSconvmat)
