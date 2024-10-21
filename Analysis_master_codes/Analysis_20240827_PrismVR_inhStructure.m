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
%%
f=26; load(fullfile(fpath{f},'PC_Result.mat'),'Result')
cd(fpath{f})
%rois={basal_ROI{f},apical_ROI{f}};
nROI=size(Result.normTraces,1);
nTau_bAP=[-20:20];
nTau={[-30:40],[-70:100],[-30:20]}; %SS, CS, dSP
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
STA_SS= STA_SS - mean(STA_SS(:,1:10),2);
F_ref=mean(STA_SS(:,-nTau_bAP(1)+[5:10]),2);
%F_ref=(tovec(imgaussfilt(Result.ref_im,1))'*tovec(Result.ftprnt)/Result.SpikeHeight_fit(1))';

SilentPeriod=ones(1,size(Result.traces,2));
sp_time=find(max(Result.spike,[],1))';
sp_na=sum((find(max(Result.spike,[],1))'+[-10:150])<0 | (find(max(Result.spike,[],1))'+[-10:150])>size(Result.traces,2),2)==0;
SilentPeriod(sp_time(sp_na)+[-10:150])=NaN;
t_fit=find(~isnan(SilentPeriod) & Result.Blue==0);

NormalizedTrace=(Result.normTraces)./F_ref;
% lwpass=NaN(size(NormalizedTrace));
% lwpass(:,t_fit)=NormalizedTrace(:,t_fit);
% lwpass=movmedian(lwpass,30000,2,'omitnan');
% NormalizedTrace=NormalizedTrace-lwpass;
NormalizedTrace_dirt=NormalizedTrace;
NormalizedTrace_dirt(Result.dirtTrace>0)=NaN;
%NormalizedTrace_dirt(:,17000:18000)=NaN;
NormalizedTrace_dirt(:,Result.motionReject)=NaN;
%NormalizedTrace_ch=cellfun(@(x) x./F_ref,Result.norm_trace_check,'UniformOutput',false);

STA_CSmat=reshape(NormalizedTrace_dirt(:,CS_s'+nTau{2}),nROI,[],length(nTau{2}));
STA_SSmat=reshape(NormalizedTrace_dirt(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1}));
STA_dSPmat=reshape(NormalizedTrace_dirt(:,dSP_s'+nTau{3}),nROI,[],length(nTau{3}));

%STA_CSmat_ch=cellfun(@(x) reshape(x(:,CS_s'+nTau{2}),nROI,[],length(nTau{2})),Result.norm_trace_check,'UniformOutput',false);
%STA_SSmat_ch=cellfun(@(x) reshape(NormalizedTrace(:,bAP_s'+nTau{1}),nROI,[],length(nTau{1})),Result.norm_trace_check,'UniformOutput',false);
%STA_dSPmat_ch=cellfun(@(x) reshape(NormalizedTrace(:,dSP_s'+nTau{3}),nROI,[],length(nTau{3})),Result.norm_trace_check,'UniformOutput',false);

%%
figure(5); clf; cmap=distinguishable_colors(6); ax1=[];
noi=setdiff([1:nROI],[BadROI{f}]);
noi_dist=ismember(Result.dist_order,noi);
%rois={basal_ROI{f},apical_ROI{f}};
cax=[-0.5 2];%*0.01;
tiledlayout(2,3);
ax3=nexttile([1 1]);
imagesc(nTau{1},[1:nROI],squeeze(mean(STA_SSmat(dist_order(noi_dist),:,:)-mean(STA_SSmat(dist_order(noi_dist),:,1:10),3,'omitnan'),2,'omitnan')),cax)
title('Simple spike')
ax4=nexttile([1 1]);
imagesc(nTau{2},[1:nROI],squeeze(mean(STA_CSmat(dist_order(noi_dist),:,:)-mean(STA_CSmat(dist_order(noi_dist),:,1:10),3,'omitnan'),2,'omitnan')),cax); %shading interp;
title('Complex spike')
ax5=nexttile([1 1]);
imagesc(nTau{3},[1:nROI],squeeze(mean(STA_dSPmat(dist_order(noi_dist),:,:)-mean(STA_dSPmat(dist_order(noi_dist),:,1:10),3,'omitnan'),2,'omitnan')),cax)
title('Dendritic spike')
colormap(turbo); linkaxes([ax3 ax4 ax5],'xy')

ax1=[ax1 nexttile([1 1])];
l=plot(nTau{1},squeeze(mean(STA_SSmat(dist_order(noi_dist),:,:)-mean(STA_SSmat(dist_order(noi_dist),:,1:10),3,'omitnan'),2,'omitnan'))); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
ax1=[ax1 nexttile([1 1])];
l=plot(nTau{2},squeeze(mean(STA_CSmat(dist_order(noi_dist),:,:)-mean(STA_CSmat(dist_order(noi_dist),:,1:10),3,'omitnan'),2,'omitnan'))); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
ax1=[ax1 nexttile([1 1])];
l=plot(nTau{3},squeeze(mean(STA_dSPmat(dist_order(noi_dist),:,:)-mean(STA_dSPmat(dist_order(noi_dist),:,1:10),3,'omitnan'),2,'omitnan'))); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
linkaxes([ax1],'xy')

figure(6); clf;
blueOn_sp=mean(Result.SpClass(:,Result.Blue>0 & Result.motionReject==0),2,'omitnan')*1000;
blueOff_sp=mean(Result.SpClass(:,Result.Blue==0 & Result.motionReject==0),2,'omitnan')*1000;
b=bar([1:3],[blueOn_sp blueOff_sp]','BarLayout','grouped','FaceColor','flat');
b(1).CData = repmat([0 0.6 1],3,1); b(2).CData = repmat([0.3 0.3 0.3],3,1);
set(gca,'XTick',[1 2 3],'XTickLabel',{'SS','CS','dSP'},'yscale','log')
ylabel('Firing rate (Hz)');
legend({'Blue On','Blue Off'})

%% compare SS waveform upon blue

f=25; bound=6;
nROI=size(Result.normTraces,1);
nTau={[-30:40],[-70:100],[-30:20]}; %SS, CS, dSP
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
cax=[-0.5 1.5];
ax1=[ax1 nexttile([1 1])];
STAss_Blueon=squeeze(mean(STA_SSmat(dist_order(noi_dist),bAP_s_on,:)-mean(STA_SSmat(dist_order(noi_dist),bAP_s_on,1:20),3,'omitnan'),2,'omitnan'));
STAss_Blueoff=squeeze(mean(STA_SSmat(dist_order(noi_dist),bAP_s_off,:)-mean(STA_SSmat(dist_order(noi_dist),bAP_s_off,1:20),3,'omitnan'),2,'omitnan'));
% STAss_Blueon=squeeze(mean(STA_SSmat(dist_order(noi_dist),bAP_s_on,:),2,'omitnan'));
% STAss_Blueoff=squeeze(mean(STA_SSmat(dist_order(noi_dist),bAP_s_off,:),2,'omitnan'));
STAss_Blueon=STAss_Blueon-(mean(STAss_Blueon(:,-nTau{1}(1)+[-10:-3]),2)-mean(STAss_Blueoff(:,-nTau{1}(1)+[-10:-3]),2));
%STAss_Blueon=STAss_Blueon-(max(STAss_Blueon,[],2)-max(STAss_Blueoff,[],2));

imagesc(nTau{1},[1:length(noi)],STAss_Blueon,cax)
title('SS during blue on'); colormap(turbo);

ax1=[ax1 nexttile([1 1])];
imagesc(nTau{1},[1:length(noi)],STAss_Blueoff,cax)
title('SS during blue off'); colormap(turbo);

ax1=[ax1 nexttile([1 1])];
imagesc(nTau{1},[1:length(noi)],movmean(STAss_Blueon-STAss_Blueoff,7,2))
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
l=plot(nTau{1},movmean(STAss_Blueon-STAss_Blueoff,7,2)); %shading interp;
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(turbo(sum(noi_dist)),2));
linkaxes([ax1],'x')
xlim([-10 40])

%% Temporal correlation
rois={basal_ROI{f},apical_ROI{f},ref_ROI{f}}; autocorrBA=[]; crosscorrBA_Blue=[];
nTau_silent=[-20:100];
all_spike=find(max(Result.SpClass([1 2 3],:),[],1));

BlueOnsetTime=find(Result.Blue>0,1);
SilentTrace=NormalizedTrace_dirt;
%SilentTrace(:,all_spike'+[-8:20])=0;
%SilentTrace(:,setdiff([1:size(Result.traces,2)],all_spike'+[-20:50]))=0;
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
title('Before Stimulation'); xlabel('Time (ms)');

ax1=[ax1 nexttile([1 1])];
plot(lagTau,crosscorrBA_Blue{2},'k'); hold all
plot(lagTau,autocorrBA_Blue{2}); hold all
legend({'Cross Basal-Apical','Auto Basal','Auto Apical','Auto Soma'})
title('During Stimulation'); xlabel('Time (ms)');
linkaxes(ax1,'x')

%%
window_sz=1000; moving_bin=1000; corrMat_moving=[];
t_slide=[1:moving_bin:size(Result.traces,2)];
for t=1:length(t_slide)-1
    t_corr=[t_slide(t):t_slide(t+1)];
    [corrMat_moving(t,:,1), t_lag]=xcorr(SilentTrace_BA(1,t_corr),SilentTrace_BA(2,t_corr),window_sz);
end

figure(8); clf;
tiledlayout(2,1);
ax1=nexttile([1 1]);
plot(t_slide(1:end-1),mean(corrMat_moving(:,find(t_lag>-4 & t_lag<4),1),2))
hold all
plot(t_slide(1:end-1),mean(corrMat_moving(:,find(t_lag>-100 & t_lag<100),1),2))
ax2=nexttile([1 1]);
imagesc(SilentTrace(Result.dist_order(noi),:),[-1 3])
linkaxes([ax1 ax2],'x')
%% Blue light unmasked subthreshold place field
nBin=150; Lap_Sub=[]; StimOn_Lap=[13:23];
all_spike=find(max(Result.SpClass([1 2],:),[],1));
Subthreshold=get_subthreshold(NormalizedTrace,all_spike,7,17);
Subthreshold_silent=Subthreshold;
Subthreshold_silent(:,all_spike'+[-8:20])=NaN;
Subthreshold_silent(:,Result.motionReject)=NaN;
Subthreshold_silent(:,Result.CStrace>0)=NaN;
Subthreshold_silent(Result.dirtTrace>0)=NaN;
bwBlue=bwlabel(Result.Blue>0);

for r=1:3 %basal, apical, soma
    Lap_Sub{r}=PlaceTrigger_average(mean(Subthreshold_silent(rois{r},:),1,'omitnan'),nBin,Result.VR,0,115,'rate');
end
Lap_FR = PlaceTrigger_average(Result.spike(1,:),nBin,Result.VR,0,115,'rate'); %total trace
Lap_CS = PlaceTrigger_average(Result.SpClass(2,:),nBin,Result.VR,0,115,'rate'); %total trace
Sp_tr=Result.spike(1,:); Sp_tr(Result.Blue==0)=NaN; %during blue on
Lap_FR_BlueOn=PlaceTrigger_average(Sp_tr,nBin,Result.VR,0,115,'rate'); %soma during blue on
Lap_blue=PlaceTrigger_average(Result.Blue,nBin,Result.VR,0,115,'rate'); %Blue

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
    imagesc(ringmovMean(Lap_Sub{r},5)); hold all
    plot(ones(1,length(StimOn_Lap))*3,StimOn_Lap,'marker','>','MarkerFaceColor',[0 0.5 1],'MarkerEdgeColor',[0 0 0])
    title([title_str{r} ' subthreshold'])
end
legend({'Stimulation Lap'},'Location','northwest')
arrayfun(@(x) colormap(x,gen_colormap([0 0.5 1; 1 1 1; 1 0 0])),ax3)

nexttile([1 1]);
Lap_sub_mean=cell2mat(cellfun(@(x) mean(ringmovMean(x(3:12,:),5),1,'omitnan'),Lap_Sub,'UniformOutput',false)');
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
%% Upward going and downward going
binEdge=[-5:0.1:5]; up_thres=1.7; dw_thres=1.5; blueTargetF=mean(reshape(NormalizedTrace(23,bAP_ref(find(ind))'+[6:10]),sum(ind),[]),[1 2]);
MedianF=zeros(2,max(bwBlue),3); %baseline during stimulation on and off
MedianVRposition=zeros(2,max(bwBlue)); %baseline during stimulation on and off
figure(12); clf; tiledlayout(3,3)
upward_pk=cell(3,2); downward_pk=cell(3,2); smooth_window=1200;
upward_trace=zeros(3,size(Result.traces,2)); downward_trace=zeros(3,size(Result.traces,2));
for r=1:3
    %before blue on
    t=[1:find(bwBlue>0,1)];
    dSub=mean(Subthreshold_silent(rois{r},t),1)-movmedian(mean(Subthreshold_silent(rois{r},t),1,'omitnan'),smooth_window,'omitnan');
    dSub=movmean(dSub,5);

    [~, upward_trace(r,t)]=detect_transient(dSub,[up_thres 1],zeros(1,length(t)));
    [~, downward_trace(r,t)]=detect_transient(-dSub,[dw_thres 1],zeros(1,length(t)));

    %during blue on

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
end

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
noi=setdiff([1:nROI],[2 3 4 5 10 28 29]); noi_dist=ismember(Result.dist_order,noi);
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
    title([title_str{r} ', ext., PC' num2str(n) ', ' num2str(D_silentOn(n)/sum(D_silentOn)*100,2) '%'])
end

for n=1:6 %Silent during blue off
    nexttile([1 1])
    imshow2(max(double(Result.ftprnt(:,:,Result.dist_order(noi_dist))>0).*reshape(rescale(V_silentOff(:,n)),1,1,[]),[],3),[])
    colormap('turbo')
    title([title_str{r} ', ext., PC' num2str(n) ', ' num2str(D_silentOff(n)/sum(D_silentOff)*100,2) '%'])
end

for n=1:6 %Pre-spike during blue on
    nexttile([1 1])
    imshow2(max(double(Result.ftprnt(:,:,Result.dist_order(noi_dist))>0).*reshape(rescale(V_preSpOn(:,n)),1,1,[]),[],3),[])
    colormap('turbo')
    title([title_str{r} ', ext., PC' num2str(n) ', ' num2str(D_preSpOn(n)/sum(D_preSpOn)*100,2) '%'])
end

for n=1:6 %Pre-spike during blue off
    nexttile([1 1])
    imshow2(max(double(Result.ftprnt(:,:,Result.dist_order(noi_dist))>0).*reshape(rescale(V_preSpOff(:,n)),1,1,[]),[],3),[])
    colormap('turbo')
    title([title_str{r} ', ext., PC' num2str(n) ', ' num2str(D_preSpOff(n)/sum(D_preSpOff)*100,2) '%'])
end

Vs=[V_silentOn(:,1:n_silentOn) V_silentOff(:,1:n_silentOff) V_preSpOn(:,1:n_preSpOn) V_preSpOff(:,1:n_preSpOff)];
Ds=[D_silentOn(1:n_silentOn)'/sum(D_silentOn) D_silentOff(1:n_silentOff)'/sum(D_silentOff) ...
    D_preSpOn(1:n_preSpOn)'/sum(D_preSpOn) D_preSpOff(1:n_preSpOff)'/sum(D_preSpOff)];
Vs=rescale2(Vs.*sign(skewness(Vs,0,1)),1);
distances = pdist(Vs', 'euclidean');  % Euclidean distance
Z = linkage(distances, 'average');
num_clusters = 7;
cluster_indices = cluster(Z, 'maxclust', num_clusters);
Vs_reduce = zeros(size(Vs, 1), num_clusters); cluster_weight=zeros(1,num_clusters);
weight_indices=cluster_indices.*Ds';
for i = 1:num_clusters
    Vs_reduce(:, i) = mean(Vs(:, find(cluster_indices == i)),2);
    %Vs_reduce(:, i) = Vs(:, find(cluster_indices == i,1));
    cluster_weight(i)=mean(weight_indices(find(cluster_indices == i)));
end

[~,sort_cluster]=sort(cluster_weight,'descend');
Vs_reduce=Vs_reduce(:,sort_cluster);
% Alternatively, you can plot a dendrogram to visualize the clustering
%imagesc(rescale2(Vs_reduce(:,sort_cluster),1))
figure(15); clf; tiledlayout(1,num_clusters);
for n=1:num_clusters %Silent during blue off
    nexttile([1 1])
    imshow2(max(double(Result.ftprnt(:,:,Result.dist_order(noi_dist))>0).*reshape(rescale(Vs_reduce(:,n)),1,1,[]),[],3),[])
    colormap('turbo')
    title([title_str{r} ', ext., PC' num2str(n) ', ' num2str(cluster_weight(sort_cluster(n))/sum(cluster_weight)*100,2) '%'])
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