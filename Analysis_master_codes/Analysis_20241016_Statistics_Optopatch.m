
clear
clc;
cd '/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/' ...
    'Prism_OptopatchData_Arrangement.xlsx'], 'Sheet1', 'B5:K175');

save_to='/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
fpath=raw(:,1);
Mouse=cell2mat(raw(:,2));
NeuronInd=cell2mat(raw(:,5));
CamType=raw(:,3);
StructureData=raw(:,10);

place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
set(0,'DefaultFigureWindowStyle','docked')


%% synchronize
isResult=zeros(1,length(fpath));
isTrace=zeros(1,length(fpath));
for f=1:length(fpath)
    f
resultfilelist=dir(fullfile(fpath{f},['OP_Result*']));
if ~isempty(resultfilelist)
    clear Result
load(fullfile(fpath{f},resultfilelist(1).name),'Result')
% save(fullfile(fpath{f},'OP_Result.mat'),'Result','-v7.3')
if contains(fpath{f}, 'BHL18TB_D2')
backupServer(fpath{f},'BHL18TB_D2','cohen_lab/Lab/Labmembers/Byung Hun Lee/Data','OP_Result.mat')
else
backupServer(fpath{f},'cohen_lab/Lab/Labmembers/Byung Hun Lee/Data','BHL18TB_D2','OP_Result.mat')    
end
isResult(f)=1;
if isfield(Result,'traces')
isTrace(f)=1;    
end
end
end

%%

isResult=zeros(1,length(fpath));
for f=85:length(fpath)
resultfilelist=dir(fullfile(fpath{f},['Result*']));
if ~isempty(resultfilelist)
load(fullfile(fpath{f},resultfilelist(1).name),'Result')
save(fullfile(fpath{f},'OP_Result.mat'),'Result')
try
backupServer(fpath{f},'BHL18TB_D2','cohen_lab/Lab/Labmembers/Byung Hun Lee/Data','OP_Result.mat')
catch
backupServer(fpath{f},'cohen_lab/Lab/Labmembers/Byung Hun Lee/Data','BHL18TB_D2','OP_Result.mat')    
end
isResult(f)=1;
end
end

%%
for f=26%:length(fpath)
load(fullfile(fpath{f},'PC_Result.mat'),'Result')    
interDendDist=[];
SkelDend = Skeletonize_dendrite(Result.ref_im,6,0.02,10);
nROI=size(Result.ftprnt,3);
for i=1:nROI
i
for j=1:nROI
[interDendDist(i,j), path]=geodesic_distance(SkelDend,get_coord(Result.ftprnt(:,:,i)),get_coord(Result.ftprnt(:,:,j)));
end
end
Result.interDendDist=interDendDist;
save(fullfile(fpath{f},'PC_Result.mat'),'Result')
end

%%
f=26;%:length(fpath)
load(fullfile(fpath{f},'PC_Result.mat'),'Result');    
figure(33); clf;
nexttile([1 1]);
imshow2(Result.ref_im,[]);
nexttile([1 1]);
plot(mean(Result.traces(ref_ROI{f},:),1))
%%

f=26; load(fullfile(fpath{f},'PC_Result.mat'),'Result')
cd(fpath{f})
rois={basal_ROI{f},apical_ROI{f},ref_ROI{f}};
nROI=size(Result.normTraces,1);
nTau_bAP=[-20:20];
nTau={[-70:50],[-70:50],[-30:20],[-70:50]}; %SS, CS, dSP, Brst
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
SS_s=[];
for s=som_spike
    isnearby=sum(ismember(s+nTau{1},som_spike))>1;
    isnearbyCS=sum(ismember(s+nTau{1},find(Result.CStrace)))>1;
    ispartCS=tr_trace(s)>0;
    if ~isnearby & ~isnan(s) & ~isnearbyCS & ~ispartCS
        SS_s=[SS_s s];
    end
end
sp_na=sum((SS_s'+nTau{1})<0 | (SS_s'+nTau{1})>size(Result.traces,2),2)==0;
SS_s=SS_s(sp_na);

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


STA_SS=squeeze(mean(reshape(Result.normTraces(:,bAP_ref'+nTau_bAP),nROI,[],length(nTau_bAP)),2));
STA_SS= STA_SS - mean(STA_SS(:,1:10),2);
F_ref=mean(STA_SS(:,-nTau_bAP(1)+[7:11]),2);
%F_ref=(tovec(imgaussfilt(Result.ref_im,1))'*tovec(Result.ftprnt)/Result.SpikeHeight_fit(1))';

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
NormalizedTrace_dirt(Result.dirtTrace>0)=NaN;
%NormalizedTrace_dirt(:,17000:18000)=NaN;
NormalizedTrace_dirt(:,Result.motionReject)=NaN;
NormalizedTrace_ch=cellfun(@(x) x./F_ref,Result.norm_trace_check,'UniformOutput',false);
NormalizedTrace_ch{1}(:,Result.motionReject)=NaN; NormalizedTrace_ch{2}(:,Result.motionReject)=NaN;
NormalizedTrace_ch{1}(Result.dirtTrace>0)=NaN; NormalizedTrace_ch{2}(Result.dirtTrace>0)=NaN;

STA_CSmat=reshape(NormalizedTrace_dirt(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));
STA_SSmat=reshape(NormalizedTrace_dirt(:,SS_s'+nTau{1}),nROI,[],length(nTau{1}));
STA_dSPmat=reshape(NormalizedTrace_dirt(:,dSP_s'+nTau{3}),nROI,[],length(nTau{3}));
STA_BSmat=reshape(NormalizedTrace_dirt(:,BS_s'+nTau{4}),nROI,[],length(nTau{4}));

Spike_BSmat=squeeze(reshape(Result.spike(1,BS_s'+nTau{4}),1,[],length(nTau{4})));
Spike_CSmat=squeeze(reshape(Result.spike(1,CS_s'+nTau{2}),1,[],length(nTau{2})));

STA_CSmat_ch=cellfun(@(x) reshape(x(:,CS_s'+nTau{2}),nROI,[],length(nTau{2})),NormalizedTrace_ch,'UniformOutput',false);
STA_SSmat_ch=cellfun(@(x) reshape(x(:,SS_s'+nTau{1}),nROI,[],length(nTau{1})),NormalizedTrace_ch,'UniformOutput',false);
STA_dSPmat_ch=cellfun(@(x) reshape(x(:,dSP_s'+nTau{3}),nROI,[],length(nTau{3})),NormalizedTrace_ch,'UniformOutput',false);
STA_BSmat_ch=cellfun(@(x) reshape(x(:,BS_s'+nTau{4}),nROI,[],length(nTau{4})),NormalizedTrace_ch,'UniformOutput',false);

SS_offind=find(Result.Blue(SS_s)==0);
CS_offind=find(Result.Blue(CS_s)==0);
BS_offind=find(Result.Blue(BS_s)==0);
%% inter distance matrix

% set the ROIs, pre-spike window
noi=setdiff([1:size(Result.traces,1)],BadROI{f});
noi_dist=ismember(Result.dist_order,noi);
Soma_ind=find(Result.dist_order(noi_dist)==1);

interDendDist=Result.interDendDist(Result.dist_order(noi_dist),Result.dist_order(noi_dist));
%%
% set the coordinates
coord_1d=dim_reduce(get_coord(Result.ftprnt));
coord_1d=coord_1d-coord_1d(1);
coord_BAsign=sign(-coord_1d(Result.dist_order(noi_dist)));
dist_from_Soma=interDendDist(:,Soma_ind)*1.17.*coord_BAsign;

% extract the pre-spike subthresholds of SS, CS, and BS
pre_spike_window=[-5:-1];
SSpresub=mean(STA_SSmat(Result.dist_order(noi_dist),SS_offind,-nTau{1}(1)+pre_spike_window),3,'omitnan');
SSdVdx=[]; for s=1:size(SSpresub,2); p = polyfit(dist_from_Soma, SSpresub(:,s), 1); SSdVdx(s)=p(1); end

BSpresub=mean(STA_BSmat(Result.dist_order(noi_dist),BS_offind,-nTau{4}(1)+pre_spike_window),3,'omitnan');
BSdVdx=[]; for s=1:size(BSpresub,2); p = polyfit(dist_from_Soma, BSpresub(:,s), 1); BSdVdx(s)=p(1); end

CSpresub=mean(STA_CSmat(Result.dist_order(noi_dist),CS_offind,-nTau{2}(1)+pre_spike_window),3,'omitnan');
CSdVdx=[]; for s=1:size(CSpresub,2); p = polyfit(dist_from_Soma, CSpresub(:,s), 1); CSdVdx(s)=p(1); end

% extract the bAP amplitude of SS, n th spike of CS, and BS
SSspAmp=max(STA_SSmat(Result.dist_order(noi_dist),SS_offind,-nTau{1}(1)+[0:4]),[],3);
SSspdelay=get_delay(squeeze(mean(STA_SSmat(Result.dist_order(noi_dist),SS_offind,-nTau{1}(1)+[-2:4]),2,'omitnan')),100);
SSspdelay=SSspdelay-SSspdelay(find(Result.dist_order(noi_dist)==1));
SSspAmp_kink=SSspAmp-SSpresub;
SSspAmp_norm=SSspAmp./SSspAmp(Soma_ind,:);

CS_spikeind=find(Spike_CSmat(CS_offind,:)');
STA_CSmat_perm=permute(STA_CSmat(:,CS_offind,:),[1 3 2]);
STA_CSmat_perm_kink=permute(STA_CSmat(Result.dist_order(noi_dist),CS_offind,:)-CSpresub,[1 3 2]);

STA_CSmat_perm_crop=reshape(STA_CSmat_perm(Result.dist_order(noi_dist),CS_spikeind'+[-2:4]'),sum(noi_dist),7,[]);
CSspAmp=max(STA_CSmat_perm_crop,[],2);
CSspAmp=squeeze(CSspAmp);
CSspAmp_kink=max(reshape(STA_CSmat_perm_kink(:,CS_spikeind'+[0:4]'),sum(noi_dist),5,[]),[],2);
CSspAmp_kink=squeeze(CSspAmp_kink);
CSspAmp_norm=CSspAmp./CSspAmp(Soma_ind,:);

BS_spikeind=find(Spike_BSmat(BS_offind,:)');
STA_BSmat_perm=permute(STA_BSmat(:,BS_offind,:),[1 3 2]);
STA_BSmat_perm_kink=permute(STA_BSmat(Result.dist_order(noi_dist),BS_offind,:)-BSpresub,[1 3 2]);

STA_BSmat_perm_crop=reshape(STA_BSmat_perm(Result.dist_order(noi_dist),BS_spikeind'+[-2:4]'),sum(noi_dist),7,[]);
BSspAmp=max(reshape(STA_BSmat_perm(Result.dist_order(noi_dist),BS_spikeind'+[0:4]'),sum(noi_dist),5,[]),[],2);
BSspAmp=squeeze(BSspAmp);
BSspAmp_kink=max(reshape(STA_BSmat_perm_kink(:,BS_spikeind'+[0:4]'),sum(noi_dist),5,[]),[],2);
BSspAmp_kink=squeeze(BSspAmp_kink);
BSspAmp_norm=BSspAmp./BSspAmp(Soma_ind,:);

% labeling n th order of spike for CS and BS
Spike_CSorderMat=NaN(size(Spike_CSmat(CS_offind,:),1),max(sum(Spike_CSmat(CS_offind,:),2)));
for c=1:size(Spike_CSmat(CS_offind,:),1)
    Spike_CSorderMat(c,1:sum(Spike_CSmat(CS_offind(c),:),2))=1;
end
[CSsp_order, CSsp_event]=ind2sub([size(Spike_CSorderMat,2) size(Spike_CSorderMat,1)],find(Spike_CSorderMat'==1));
Spike_CSorderMat_ind=Spike_CSorderMat';
Spike_CSorderMat_ind(~isnan(Spike_CSorderMat_ind))=[1:sum(~isnan(Spike_CSorderMat_ind(:)))];
Spike_CSorderMat_ind=Spike_CSorderMat_ind';

Spike_BSorderMat=NaN(size(Spike_BSmat(BS_offind,:),1),max(sum(Spike_BSmat(BS_offind,:),2)));
for c=1:size(Spike_BSmat(BS_offind,:),1)
    Spike_BSorderMat(c,1:sum(Spike_BSmat(BS_offind(c),:),2))=1;
end
[BSsp_order, ~]=ind2sub([size(Spike_BSorderMat,2) size(Spike_BSorderMat,1)],find(Spike_BSorderMat'==1));
Spike_BSorderMat_ind=Spike_BSorderMat';
Spike_BSorderMat_ind(~isnan(Spike_BSorderMat_ind))=[1:sum(~isnan(Spike_BSorderMat_ind(:)))];
Spike_BSorderMat_ind=Spike_BSorderMat_ind';

% mean of n th order of spike for CS
delay_mat=NaN(sum(noi_dist),4,3);

figure(102); clf;
tiledlayout(3,4);
nexttile([1 1]);
imagesc(squeeze(mean(STA_SSmat(Result.dist_order(noi_dist),SS_offind,-nTau{1}(1)+1+[-2:4]),2,'omitnan')),[-2 10])
xlabel('Time (ms)'); ylabel('Basal to apical ROI');
title('Average of SS')
delay_mat(:,1,1)=get_delay(squeeze(mean(STA_SSmat(Result.dist_order(noi_dist),SS_offind,-nTau{1}(1)+1+[-2:4]),2,'omitnan')),100);
delay_mat(:,1,1)=delay_mat(:,1,1)-min(delay_mat(:,1,1));

for s=1:3
nexttile([1 1]);
imagesc(mean(STA_BSmat_perm_crop(:,:,BSsp_order==s),3),[-2 10])
xlabel('Time (ms)'); ylabel('Basal to apical ROI');
title(['Average of ' counting_string(s) ' spike of BS'])
delay_mat(:,s,3)=get_delay(mean(STA_BSmat_perm_crop(:,:,BSsp_order==s),3),100);
delay_mat(:,s,3)=delay_mat(:,s,3)-min(delay_mat(:,s,3));
end

for s=1:4
nexttile([1 1]);
imagesc(mean(STA_CSmat_perm_crop(:,:,CSsp_order==s),3),[-2 10])
xlabel('Time (ms)'); ylabel('Basal to apical ROI');
title(['Average of ' counting_string(s) ' spike of CS'])
delay_mat(:,s,2)=get_delay(mean(STA_CSmat_perm_crop(:,:,CSsp_order==s),3),100);
delay_mat(:,s,2)=delay_mat(:,s,2)-min(delay_mat(:,s,2));
end
colormap(turbo)

nexttile([1 3]); cmap_SS=[0 0 0]; cmap_CS=[1 0 0; 1 0.2 0.2; 1 0.4 0.4; 1 0.6 0.6];
cmap_BS=[0 0 1; 0.2 0.2 1; 0.4 0.4 1];
[~, showind]=sort(dist_from_Soma,'ascend');
plot(dist_from_Soma(showind),delay_mat(showind,1,1),'.-','color',cmap_SS); hold all
for s=1:4; plot(dist_from_Soma(showind),delay_mat(showind,s,2),'.-','color',cmap_CS(s,:)); end
for s=1:3; plot(dist_from_Soma(showind),delay_mat(showind,s,3),'.-','color',cmap_BS(s,:)); end
legend({'SS','CS 1st','CS 2nd','CS 3rd','CS 4th','BS 1st','BS 2nd','BS 3rd'})
xlabel('Distance from soma (\mum)'); ylabel('bAP delay (ms)')

conduction_Speed=NaN(size(delay_mat,2),size(delay_mat,3));
for s1=1:size(delay_mat,2)
for s2=1:size(delay_mat,3)
p = polyfit(abs(dist_from_Soma), delay_mat(:,s1,s2), 1);
conduction_Speed(s1,s2)=1/p(1);
end
end
nexttile([1 1])
scatter([1:sum(~isnan(conduction_Speed(:)))],conduction_Speed(~isnan(conduction_Speed)),15,[cmap_SS; cmap_CS; cmap_BS],'filled')
set(gca,'XTick',[1:sum(~isnan(conduction_Speed(:)))],'XTickLabel',{'SS','CS 1st','CS 2nd','CS 3rd','CS 4th','BS 1st','BS 2nd','BS 3rd'})
ylabel('Conduction velocity (\mu/ms)')

% measure the slope between pre-spike subthreshold and bAP amplitude
SS_AmpsubCorr=[]; CS_AmpsubCorr=[]; BS_AmpsubCorr=[];
for n=1:sum(noi_dist);
    rmvnan=find(~isnan(SSpresub(n,:)) & ~isnan(SSspAmp(n,:)));
    SS_AmpsubCorr(n,1) = corr(SSpresub(n,rmvnan)', SSspAmp_kink(n,rmvnan)');
    for s=1:4
        sp_ind=find(Spike_CSorderMat(:,s)==1);
        Ampind=Spike_CSorderMat_ind(sp_ind,s);
        sp_ind2=find(~isnan(CSpresub(n,sp_ind)) & ~isnan(CSspAmp(n,Ampind)));
        CS_AmpsubCorr(n,s) = corr(CSpresub(n,sp_ind(sp_ind2))', CSspAmp_kink(n,Ampind(sp_ind2))'); 
    end
    for s=1:3
        sp_ind=find(Spike_BSorderMat(:,s)==1);
        Ampind=Spike_BSorderMat_ind(sp_ind,s);
        sp_ind2=find(~isnan(BSpresub(n,sp_ind)) & ~isnan(BSspAmp(n,Ampind)));
        BS_AmpsubCorr(n,s) = corr(BSpresub(n,sp_ind(sp_ind2))', BSspAmp_kink(n,Ampind(sp_ind2))');
    end
end

% measure the delay

figure(100); clf; cmap=turbo(sum(noi_dist)); ax1=[];

for n=1:size(SSspAmp,1)
    ax1=[ax1 nexttile([1 1])];
plot(SSpresub(n,:)',SSspAmp(n,:)','.','Color',cmap(n,:));
xlabel('Pre-spike subthreshold (\DeltaF/F_{ref})')
ylabel('Spike amplitude (\DeltaF/F_{ref})')
title(['ROI #' num2str(n) ', Distance from soma: ' num2str(interDendDist(n,Soma_ind)*1.17,3) '\mum'])
end
linkaxes(ax1,'xy');
axis tight

figure(101); clf; cmap=distinguishable_colors(6);
nexttile([1 1])
[~, showind]=sort(dist_from_Soma,'ascend');
h(1)=errorbar_shade(dist_from_Soma(showind)',mean(SSspAmp_norm(showind,:),2,'omitnan'),std(SSspAmp_norm(showind,:),0,2,'omitnan'),cmap(1,:));
hold all
sp_ind=find(CSsp_order==1);
h(2)=errorbar_shade(dist_from_Soma(showind)',mean(CSspAmp_norm(:,sp_ind),2,'omitnan'),std(CSspAmp_norm(:,sp_ind),0,2,'omitnan'),cmap(2,:));
sp_ind=find(BSsp_order==1);
h(3)=errorbar_shade(dist_from_Soma(showind)',mean(BSspAmp_norm(:,sp_ind),2,'omitnan'),std(BSspAmp_norm(:,sp_ind),0,2,'omitnan'),cmap(3,:));
legend(h,{'SS','CS','BS'})
xlabel('Distance from soma (\mum)')
ylabel('Normalized bAP amplitude')
ylim([0.3 1.2]); xlim([min(dist_from_Soma) max(dist_from_Soma)]);

nexttile([1 1])
cmap=distinguishable_colors(4); h=[];
for s=1:4
sp_ind=find(CSsp_order==s);
%h(s)=errorbar_shade(dist_from_Soma(showind)',mean(CSspAmp_norm(:,sp_ind),2,'omitnan'),std(CSspAmp_norm(:,sp_ind),0,2,'omitnan'),cmap(s,:)); hold all
h(s)=plot(dist_from_Soma(showind)',mean(CSspAmp_norm(:,sp_ind),2,'omitnan'),'color',cmap(s,:)); hold all
end
legend(h,arrayfun(@num2str, [1:4], 'UniformOutput', false))
title(['Spike attenuation of complex spike'])
xlabel('Distance from soma (\mum)')
ylabel('Normalized bAP amplitude')
ylim([0.3 1.2]); xlim([min(dist_from_Soma) max(dist_from_Soma)]);

nexttile([1 1])
h=[];
for s=1:3
sp_ind=find(BSsp_order==s);
%h(s)=errorbar_shade(dist_from_Soma(showind)',mean(CSspAmp_norm(:,sp_ind),2,'omitnan'),std(CSspAmp_norm(:,sp_ind),0,2,'omitnan'),cmap(s,:)); hold all
h(s)=plot(dist_from_Soma(showind)',mean(BSspAmp_norm(:,sp_ind),2,'omitnan'),'color',cmap(s,:)); hold all
end
legend(h,arrayfun(@num2str, [1:3], 'UniformOutput', false))
title(['Spike attenuation of burst spike'])
xlabel('Distance from soma (\mum)')
ylabel('Normalized bAP amplitude')
ylim([0.3 1.2]); xlim([min(dist_from_Soma) max(dist_from_Soma)]);

nexttile([1 1]); cmap=distinguishable_colors(6);
h(1)=errorbar_shade(dist_from_Soma(showind)',mean(SSspAmp_kink(showind,:),2,'omitnan'),std(SSspAmp_kink(showind,:),0,2,'omitnan'),cmap(1,:));
hold all
sp_ind=find(CSsp_order==1);
h(2)=errorbar_shade(dist_from_Soma(showind)',mean(CSspAmp_kink(:,sp_ind),2,'omitnan'),std(CSspAmp_kink(:,sp_ind),0,2,'omitnan'),cmap(2,:));
sp_ind=find(BSsp_order==1);
h(3)=errorbar_shade(dist_from_Soma(showind)',mean(BSspAmp_kink(:,sp_ind),2,'omitnan'),std(BSspAmp_kink(:,sp_ind),0,2,'omitnan'),cmap(3,:));
legend(h,{'SS','CS','BS'})
xlabel('Distance from soma (\mum)')
ylabel('bAP amplitude (kink)')
xlim([min(dist_from_Soma) max(dist_from_Soma)]);

nexttile([1 1])
cmap=distinguishable_colors(4); h=[];
for s=1:4
sp_ind=find(CSsp_order==s);
%h(s)=errorbar_shade(dist_from_Soma(showind)',mean(CSspAmp_norm(:,sp_ind),2,'omitnan'),std(CSspAmp_norm(:,sp_ind),0,2,'omitnan'),cmap(s,:)); hold all
h(s)=plot(dist_from_Soma(showind)',mean(CSspAmp_kink(:,sp_ind),2,'omitnan'),'color',cmap(s,:)); hold all
end
legend(h,counting_string([1:4]))
title(['Spike attenuation of complex spike'])
xlabel('Distance from soma (\mum)')
ylabel('bAP amplitude (kink)')
xlim([min(dist_from_Soma) max(dist_from_Soma)]);

nexttile([1 1])
h=[];
for s=1:3
sp_ind=find(BSsp_order==s);
%h(s)=errorbar_shade(dist_from_Soma(showind)',mean(CSspAmp_norm(:,sp_ind),2,'omitnan'),std(CSspAmp_norm(:,sp_ind),0,2,'omitnan'),cmap(s,:)); hold all
h(s)=plot(dist_from_Soma(showind)',mean(BSspAmp_kink(:,sp_ind),2,'omitnan'),'color',cmap(s,:)); hold all
end
legend(h,counting_string([1:3]))
title(['Spike attenuation of burst spike'])
xlabel('Distance from soma (\mum)')
ylabel('bAP amplitude (kink)')
xlim([min(dist_from_Soma) max(dist_from_Soma)]);

nexttile([1 1])
h=[];
plot(dist_from_Soma(showind),SS_AmpsubCorr(showind,:)); hold all
plot(dist_from_Soma(showind),CS_AmpsubCorr(showind,1))
plot(dist_from_Soma(showind),BS_AmpsubCorr(showind,1))
xlabel('Distance from soma (\mum)')
ylabel(['Correlation between' newline 'Pre-spike sub and Spike Amp. (kink)'])
xlim([min(dist_from_Soma) max(dist_from_Soma)]);
legend(h,{'SS','CS','BS'})

nexttile([1 1])
h=[];
plot(dist_from_Soma(showind),CS_AmpsubCorr(showind,:))
xlabel('Distance from soma (\mum)')
ylabel(['Correlation between' newline 'Pre-spike sub and Spike Amp. (kink)'])
xlim([min(dist_from_Soma) max(dist_from_Soma)]);
legend(counting_string([1:4]))
title('Complex spike')

nexttile([1 1])
h=[];
plot(dist_from_Soma(showind),BS_AmpsubCorr(showind,:))
xlabel('Distance from soma (\mum)')
ylabel(['Correlation between' newline 'Pre-spike sub and Spike Amp. (kink)'])
xlim([min(dist_from_Soma) max(dist_from_Soma)]);
legend(counting_string([1:3]))
title('Burst spike')

%%


