clear; clc;
[~, ~, raw] = xlsread(['/Volumes/cohen_lab/Lab/Labmembers' ...
    '/Byung Hun Lee/Data/PrismPCdata_Arrangement.xlsx'], 'Sheet1', 'C5:Z31');
fpath=raw(:,1);
foi=[1 4 5 6 8 10 11 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27];

ref_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,9),'UniformOutput',false);

oblique_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,10),'UniformOutput',false);
PeriSoma_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,11),'UniformOutput',false);
basal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,12),'UniformOutput',false);
apical_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,13),'UniformOutput',false);
distal_ROI=cellfun(@(x) (str2num(num2str(x))),raw(:,14),'UniformOutput',false);

fpath=raw(:,1)';
StructureData=raw(:,8);
BadROI=cellfun(@(x) (str2num(num2str(x))),raw(:,17),'UniformOutput',false);
EndFrame=cell2mat(raw(:,15));
ifmotionReject=cell2mat(raw(:,16));
ifdirtRemov=cell2mat(raw(:,18));
Pixelsize=cell2mat(raw(:,6));
save_figto='/Volumes/BHL_WD18TB/PP72_PlaceCellResults';
place_bin=150; time_segment=15000; overlap=200;
alignedMovFN = {'STA_Mat_SS','STA_Mat_CS','STA_Mat_dSP'};
bound=6;
title_str={'Basal','Apical','Peri-Soma'};
PlaceFieldList=cellfun(@(x) (str2num(num2str(x))),raw(:,21),'UniformOutput',false);
PlaceFieldBin=cellfun(@(x) (str2num(num2str(x))),raw(:,22),'UniformOutput',false);
set(0,'DefaultFigureWindowStyle','docked')
%%
f=26;
load(fullfile(fpath{f},'PC_Result.mat'))
   f
load(fullfile(fpath{f},'PC_Result.mat'))
Corrmat=[]; normTrace=[];
normTr=Result.normTraces./Result.F0_PCA;
normTr=normTr(Result.dist_order,:);
SubTr=get_subthreshold(normTr,Result.spike(1,:),7,17);

normTrace{1}=Result.norm_trace_check{1}./Result.F0_PCA;
normTrace{2}=Result.norm_trace_check{2}./Result.F0_PCA;
normTrace=cellfun(@(x) x(Result.dist_order,:),normTrace,'UniformOutput',false);
SubTrace=cellfun(@(x) get_subthreshold(x,Result.spike(1,:),7,17),normTrace,'UniformOutput',false);

if ifmotionReject(f)
SubTr(:,Result.motionReject)=NaN;
SubTrace{1}(:,Result.motionReject)=NaN;
SubTrace{2}(:,Result.motionReject)=NaN;
end

if ifdirtRemov(f)
    dirtTrace=Result.dirtTrace(Result.dist_order,:);
SubTr(dirtTrace>0)=NaN;
SubTrace{1}(dirtTrace>0)=NaN;
SubTrace{2}(dirtTrace>0)=NaN;
end

branchlabel=Result.BranchLabel;
branchlabel=branchlabel(Result.dist_order);
for b=1:max(branchlabel)
bind=find(branchlabel==b);
Corrmat{b}=get_corrMat(SubTrace{1}(bind,:),SubTrace{2}(bind,:));
end

SubTrace_filter=pcafilterTrace(SubTrace,5);
LapSub=PlaceTrigger_average(SubTrace_filter,150,Result.VR,-0.002,115); %total trace
LapFR=PlaceTrigger_average(Result.spike(1,:),150,Result.VR,-0.002,115)*1000; %total trace

% l=repmat(ringmovMean(LapSub,5),1,3);
% imshow_patch(l(:,120+[-75:75],:))
% colormap(turbo)

StimLap=unique(Result.VR(8,:).*(Result.Blue>0));
Show_t=[find(Result.VR(8,:)==StimLap(end)+1,1):find(Result.VR(8,:)==StimLap(end)+1,1)+140000];
show_lap=unique(Result.VR(8,Show_t));
%%
figure(11); clf; tiledlayout(3,4);
cmap=[0. 0.1 0.5; 0.5 0.1 0]; cmap_line=[0.2 0.6 1; 1 0.6 0.2]; h=[];
nexttile([2 4])
plot([Show_t-Show_t(1)]/1000,normTr(11,Show_t),'color',cmap_line(1,:)); hold all
plot([Show_t-Show_t(1)]/1000,SubTr(11,Show_t),'color',cmap(1,:))
sp=Result.spike(1,Show_t);
%plot(find(sp)/1000,normTr(11,Show_t(find(sp))),'r.')
plot(find(sp)/1000,0.012,'r.')

plot([Show_t-Show_t(1)]/1000,normTr(19,Show_t)-0.02,'color',cmap_line(2,:)); hold all
plot([Show_t-Show_t(1)]/1000,SubTr(19,Show_t)-0.02,'color',cmap(2,:))

plot([Show_t-Show_t(1)]/1000,Result.VR(5,Show_t)/10000-0.04,'color',[0.2 0.7 0])
axis off

PFcenter=75; nVRbin=[-50:150];
nexttile([1 1])
L=repmat(ringmovMean(LapFR(show_lap,:),7),1,3);
imagesc((PFcenter-150+nVRbin)/150*200,[1:length(show_lap)],L(:,PFcenter+150+nVRbin)); 
cb=colorbar;
cb.Label.String='Firing rate (Hz)';
xlabel('VR position (cm)');

nexttile([1 1])
L=repmat(ringmovMean(LapSub(show_lap,:,11),7),1,3);
imagesc((PFcenter-150+nVRbin)/150*200,[1:length(show_lap)],L(:,PFcenter+150+nVRbin),[-0.002 0.005])
xlabel('VR position (cm)');
cb=colorbar;
cb.Label.String='\DeltaF/F';

nexttile([1 1])
L=repmat(ringmovMean(LapSub(show_lap,:,19),7),1,3);
imagesc((PFcenter-150+nVRbin)/150*200,[1:length(show_lap)],L(:,PFcenter+150+nVRbin),[-0.002 0.005])
colormap(turbo);
xlabel('VR position (cm)');
cb=colorbar;
cb.Label.String='\DeltaF/F';

nexttile([1 1])
L=repmat(ringmovMean(cat(3,LapFR(show_lap,:),LapSub(show_lap,:,11),LapSub(show_lap,:,19)),5),1,3);

plot((PFcenter-150+nVRbin)/150*200,mean(L(:,PFcenter+150+nVRbin,1),1,'omitnan'),'k'); hold all
ylabel('Firing rate (Hz)');
yyaxis right
plot((PFcenter-150+nVRbin)/150*200,mean(L(:,PFcenter+150+nVRbin,2),1,'omitnan'),'color',cmap_line(1,:),'linestyle','-')
plot((PFcenter-150+nVRbin)/150*200,mean(L(:,PFcenter+150+nVRbin,3),1,'omitnan'),'color',cmap_line(2,:),'linestyle','-')
set(gca,'YColor','k')
ylabel('\DeltaF/F');
xlabel('VR position (cm)');

%%


%%
for f=26
    f
load(fullfile(fpath{f},'PC_Result.mat'))
Corrmat=[]; normTrace=[];
normTr=Result.normTraces./Result.F0_PCA;
normTr=normTr(Result.dist_order,:);
SubTr=get_subthreshold(normTr,Result.spike(1,:),7,17);

normTrace{1}=Result.norm_trace_check{1}./Result.F0_PCA;
normTrace{2}=Result.norm_trace_check{2}./Result.F0_PCA;
normTrace=cellfun(@(x) x(Result.dist_order,:),normTrace,'UniformOutput',false);
SubTrace=cellfun(@(x) get_subthreshold(x,Result.spike(1,:),7,17),normTrace,'UniformOutput',false);

if ifmotionReject(f)
SubTr(:,Result.motionReject)=NaN;
SubTrace{1}(:,Result.motionReject)=NaN;
SubTrace{2}(:,Result.motionReject)=NaN;
end

if ifdirtRemov(f)
    dirtTrace=Result.dirtTrace(Result.dist_order,:);
SubTr(dirtTrace>0)=NaN;
SubTrace{1}(dirtTrace>0)=NaN;
SubTrace{2}(dirtTrace>0)=NaN;
end

branchlabel=Result.BranchLabel;
branchlabel=branchlabel(Result.dist_order);
for b=1:max(branchlabel)
bind=find(branchlabel==b);
Corrmat{b}=get_corrMat(SubTrace{1}(bind,:),SubTrace{2}(bind,:));
end

figure(f); clf;
imshow_patch(Corrmat,[0.5 1])
colormap(turbo)
drawnow;
end

%% Load movie

frame_interest=539750;
mov_res_cat=[];
for j=[36 37]
load(fullfile(fpath{f},'PC_Result.mat')) % load the result file
load(fullfile(fpath{f},"output_data.mat")) 
sz=double(Device_Data{1, 3}.ROI([2 4]));
mov_mc=double(readBinMov([fpath{f} '/mc_ShutterReg' num2str(j,'%02d') '.bin'],sz(2),sz(1))); %load movie
load(fullfile(fpath{f},['mcTrace' num2str(j,'%02d') '.mat'])); % load motion correction traces

mov_res=mov_mc-mean(mov_mc,3);
mov_res = SeeResiduals(mov_res,mcTrace.xymean);
mov_res = SeeResiduals(mov_res,mcTrace.xymean.^2);
mov_res = SeeResiduals(mov_res,mcTrace.xymean(:,1).*mcTrace.xymean(:,end));
mov_res_cat=cat(3,mov_res_cat,mov_res);
end
%mov_res = mov_res.*(max(Result.bvMask,[],3)==0); % mask out blood vessels
mov_res_cat=mov_res_cat(:,:,15000+[-8000:8000]);
F0=imgaussfilt(mean(mov_mc,3)-100,5);
mov_res_dff=mov_res_cat./F0;

DendBound=bwboundaries(bwlabel(max(Result.ftprnt>0,[],3)>0));

figure(12); clf; tiledlayout(1,3);
ax1=nexttile([1 1]);
imshow2(imrotate(Result.ref_im(6:end-6,6:end-6),90)-100,[400 3000]); hold all
colormap(ax1,gray)
for b=1:size(DendBound,1)
plot(DendBound{b}(:,1)-6,sz(1)-DendBound{b}(:,2)-6,'w')
end
title('F0')
cb=colorbar;
cb.Label.String='F';

ax2=nexttile([1 1]);
bAP_df=max(-imgaussfilt(mov_res_cat(:,:,10365:10368),1),[],3);
bAP_df=bAP_df.*(max(Result.bvMask,[],3)==0);%.*(max(Result.ftprnt>0,[],3));
imshow2(imrotate(bAP_df(6:end-6,6:end-6),90),[30 120]); hold all
colormap(ax2,hot)
for b=1:size(DendBound,1)
plot(DendBound{b}(:,1)-6,sz(1)-DendBound{b}(:,2)-6,'w')
end
title('bAP')
cb=colorbar;
cb.Label.String='\DeltaF';

ax3=nexttile([1 1]);
bAP_df=mean(-imgaussfilt(mov_res_cat(:,:,7649:7722),1),3);
bAP_df=bAP_df.*(max(Result.bvMask,[],3)==0);%.*(max(Result.ftprnt>0,[],3));
imshow2(imrotate(bAP_df(6:end-6,6:end-6),90),[5 50]); hold all
colormap(ax3,hot)
for b=1:size(DendBound,1)
plot(DendBound{b}(:,1)-6,sz(1)-DendBound{b}(:,2)-6,'w')
end
title('Local subthreshold')
cb=colorbar;
cb.Label.String='\DeltaF';
linkaxes([ax1 ax2 ax3],'xy')

