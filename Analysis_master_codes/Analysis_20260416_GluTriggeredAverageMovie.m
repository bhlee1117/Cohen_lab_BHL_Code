%% make_GluVolt_movie
% Generates a 3-panel side-by-side movie:

clear; clc;
cd '/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
[~, ~, raw] = xlsread( ...
    '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement.xlsx', ...
    'Sheet1', 'C5:AA31');
load('/Volumes/BHL18TB_D2/20260203_SD_V2+iGluSNFR4/25X_transformationMatrix.mat');
fpath           = raw(:,1)';
V2moviemaxTime  = 15000;
GlumoviemaxTime = 5000;
StructureData=raw(:,6);
%% ===================== USER SETTINGS ====================================
f            = 3;        % session index into fpath{}
framesPerSeg  = 15000;   % frames per mc0x.bin file  (15 s × 1000 Hz)
Fs_volt       = 1000;    % voltage camera frame rate (Hz)

% Spatial smoothing (Gaussian sigma in pixels, per frame)
volt_smooth_sigma = 1.0;
glu_smooth_sigma  = 1.5;

% Optional: structure mask (H_v x W_v, values 0-1). Set [] to skip.
strImg_bin_volt = [];

%% ===================== LOAD DEVICE DATA =================================
load(fullfile(fpath{f}, 'output_data.mat'));
GluResult=importdata(fullfile(fpath{f},"Glu_Result.mat"));
VoltResult=importdata(fullfile(fpath{f},"Volt_Result.mat"));

% Voltage camera geometry
nCol1 = double(Device_Data{3}.ROI(2));
nRow1 = double(Device_Data{3}.ROI(4));
exposuretime1 = 1 / Fs_volt;   % s per frame (= 0.001 s)

% Glutamate camera geometry
nCol2 = double(Device_Data{4}.ROI(2));
nRow2 = double(Device_Data{4}.ROI(4));
exposuretime2 = Device_Data{4}.exposuretime;   % s per frame (slower camera)

%% ===================== BUILD GLU TIME AXIS ==============================
% Rebuild cam2_trig to get glutamate frame times in seconds
cam1_vsyn = Device_Data{1,2}.Counter_Inputs(1,1).data;
start_idx = find(cam1_vsyn == max(cam1_vsyn), 1);
end_idx   = length(Device_Data{1,2}.buffered_tasks(1,2).channels(1,1).data);
seg_size  = 103 * floor(exposuretime1 / 0.001);
last_val  = cam1_vsyn(start_idx - 1);
n_to_add  = end_idx - start_idx + 1;
n_seg     = ceil(n_to_add / seg_size);
added     = repelem((last_val+1 : last_val+n_seg)', seg_size);
cam1_vsyn(start_idx:end_idx) = added(1:n_to_add);

cam2_vsyn      = Device_Data{1,2}.Counter_Inputs(1,2).data;
cam2_vsyn_trig = find(diff(cam2_vsyn) == 1) + 1;
seg_size2      = cam2_vsyn_trig(10) - cam2_vsyn_trig(9);
start_idx2     = find(cam2_vsyn == max(cam2_vsyn), 1);
end_idx2       = length(Device_Data{1,2}.buffered_tasks(1,2).channels(1,2).data);
last_val2      = cam2_vsyn(start_idx2 - 1);
n_to_add2      = end_idx2 - start_idx2 + 1;
n_seg2         = ceil(n_to_add2 / seg_size2);
added2         = repelem((last_val2+1 : last_val2+n_seg2)', seg_size2);
cam2_vsyn(start_idx2:end_idx2) = added2(1:n_to_add2);

mod488 = 200001;  mod607 = 200001;
cam2_trig = find(diff(cam2_vsyn) == 1) + 1;
cam2_trig = cam1_vsyn(cam2_trig);
cam2_trig = cam2_trig(cam2_vsyn(mod488)+2 : cam2_vsyn(end))' - (cam1_vsyn(mod607)+2);

% Glu frame times in seconds (relative to volt frame 1)
t_glu = GluResult.t_ax(1:size(GluResult.dFF_glu,2));   % cam2_trig is in volt-frame units

% Rebuild GluemovTimesegments for loading
if cam2_vsyn(end) < GlumoviemaxTime
    GluemovTimesegments = [cam2_vsyn(mod488)+2,  cam2_vsyn(end)];
else
    GluemovTimesegments = [(cam2_vsyn(mod488)+2) : GlumoviemaxTime : cam2_vsyn(end)];
    GluemovTimesegments(end+1) = cam2_vsyn(end);
end

%% ===================== PRELOAD GLU MOVIE ================================
% Glu dF/F — per-ROI masked, median-filtered, and normalized
n_kernel = 1.5;   % median filter kernel size in pixels (set as needed)
dil_rad  = 2;   % dilation radius in pixels

fprintf('Loading glutamate movie segments...\n');
mov_glu =readBinMov_BHL_multiple(fpath{f},4,1:length(t_glu),5000,'mc2');
szGlu= size(mov_glu);

mov_glu2 =SeeResiduals(mov_glu  ,GluResult.mc);
mov_glu2 =SeeResiduals(mov_glu2 ,GluResult.mc.^2);
mov_glu2 =SeeResiduals(mov_glu2 ,GluResult.mc(:,1).*GluResult.mc(:,end));

[V, D, u_sub] = get_eigvector(tovec(imgaussfilt(mov_glu2,glu_smooth_sigma))',40);
u_sub=reshape(u_sub,szGlu(1),szGlu(2),[]);
figure(2); clf; imshow2_patch(u_sub); drawnow;
n=input("PCs to regress out\n");
mov_glu2 =SeeResiduals(mov_glu2 ,V(:,n));

N_roi = size(GluResult.S_glu, 3);

%-- Normalize by median of raw mov_glu within this ROI
% Compute per-pixel mean over time within the (undilated) ROI
F0_roi=tovec(mean(mov_glu,3))'* tovec(GluResult.S_glu)./squeeze(sum(tovec(GluResult.S_glu),[1 2]));

roi_mask     = max(GluResult.S_glu> 0,[],3);               % X x Y binary
roi_mask_dil = imdilate(roi_mask, strel('disk', dil_rad));% dilated mask
mov_roi=maskBinary(mov_glu2,roi_mask_dil==0,NaN);
mov_roi= imgaussfiltnan(mov_roi, n_kernel);
F0_roi_img=zeros(size(mov_glu2(:,:,1)), 'double');
for i=1:N_roi
    roi_mask     = GluResult.S_glu(:,:,i) > 0;               % X x Y binary
    roi_mask     = imdilate(roi_mask, strel('disk', dil_rad));% dilated mask
    F0_roi_img_tmp=maskBinary(zeros(size(mov_glu2(:,:,1)), 'double'),roi_mask,F0_roi(i));
    F0_roi_img=max(cat(3,F0_roi_img,F0_roi_img_tmp),[],3);
end
dFF_roi = (mov_roi) ./ F0_roi_img .*max(GluResult.S_glu> 0,[],3);
clear mov_roi mov_glu u_sub

%% Dendritic coordinates
StructureStack=read_tiff(StructureData{f});
sz_struct=size(StructureStack(:,:,1));
[H_out, W_out] = size(GluResult.AvgGluImg);

swcfiles=dir(fullfile(fpath{f}, 'Tracing*.swc'));
if isempty(swcfiles)
swcfiles=dir(fullfile(fileparts(StructureData{f}), 'Tracing*.swc'));
end
dendritePoints=[];
for si=1:length(swcfiles)
swcname=fullfile(swcfiles(si).folder, swcfiles(si).name);
tree = load_tree(swcname);
dendritePoints{si}=[tree.X, tree.Y];
end
DendriteMask=point2img(cell2mat(dendritePoints'),3.5,sz_struct);
DendriteMask_N1_str=point2img(cell2mat(dendritePoints(1)),3,sz_struct);

DendriteMask_N1= imwarp(DendriteMask_N1_str, (VoltResult.tform_Str2Volt), 'OutputView', imref2d(size(VoltResult.ref_im)));

VoltRefIm= transformCamera_O2B(Device_Data, tformReg, VoltResult.ref_im, GluResult.AvgGluImg);
VoltRefIm_rgb = grs2rgb(VoltRefIm, gray(256), 20, 200);

%% Glutamate spike & position
GluSpike=find_spike_bh(zscore(GluResult.dFF_glu,0,2),3,1);
SumTglu=sum(GluSpike,1);
SumGlu=sum(GluSpike,2);
[~, synsort]=sort(SumGlu,'ascend');
%show_traces_spikes(GluResult.dFF_glu(synsort,:),GluSpike(synsort,:),SumTglu)

F0_volt=tovec(VoltResult.ref_im)'*tovec(VoltResult.ftprnt);
VoltTrace_Norm=VoltResult.normTraces./F0_volt';

gspTA_voltage=[]; gspMat_voltage=[];
for si=1:size(GluSpike,1)
gsp_time=GluResult.t_ax(find(GluSpike(si,:)>0));
[gsptime_volt, dist] = match_nearest(gsp_time, VoltResult.t_ax);
[gspTA_voltage(:,:,si), gspMat_voltage]=get_STA(VoltTrace_Norm,gsptime_volt,100,100);
end
gspTA_voltage=gspTA_voltage-median(gspTA_voltage,2);
Glu_coord=get_coord(GluResult.S_glu);

%% Show representative synaptic event case.
VoltLag=1000; %ms
Syn2show=7;
caxis=[-0.005 0.01];
GluTimeGap=min(diff(GluResult.t_ax));
GluLag=ceil(VoltLag/(GluTimeGap*1000));
gsp_time=GluResult.t_ax(find(GluSpike(Syn2show,:)>0));
[gsptime_volt, dist] = match_nearest(gsp_time, VoltResult.t_ax);
[~, gspMat]=get_STA(GluSpike,GluSpike(Syn2show,:)>0,GluLag,GluLag);
[~, gspMat_voltage, validglu]=get_STA(VoltTrace_Norm,gsptime_volt,VoltLag,VoltLag);
validglu=ismember(gsptime_volt,validglu);
gspMat=gspMat(:,validglu,:);
gspMat_voltage=gspMat_voltage-median(gspMat_voltage,3);
AlignedVolt_ftprnt = transformCamera_O2B(Device_Data, tformReg, VoltResult.ftprnt, GluResult.AvgGluImg);

nROI=size(VoltResult.ftprnt,3);
SkelDend = DendriteMask_N1;
interDendDist=[];
% figure(99); clf;
for i=1:nROI
        [interDendDist(i), p]=geodesic_distance(SkelDend,get_coord(AlignedVolt_ftprnt(:,:,i)),get_coord(GluResult.S_glu(:,:,Syn2show)));
        % nexttile([1 1]);
        % imshow2(im_merge(cat(3,mat2gray(max(AlignedVolt_ftprnt,[],3)),p),[1 1 1; 1 0 0]),[])
end
[InterGludist_sorted, voltROIsorted]=sort(interDendDist,'ascend');
GluFtprntDistFromSyn=[];
for si=[1:size(GluResult.S_glu,3)]
[GluFtprntDistFromSyn(si) p2]=geodesic_distance(SkelDend,Glu_coord(si,:),Glu_coord(Syn2show,:));
%nexttile([1 1]); imshow2(im_merge(cat(3,mat2gray(max(AlignedVolt_ftprnt,[],3)),p2),[1 1 1; 1 0 0]),[])
end
[~, GluROIsorted]=sort(GluFtprntDistFromSyn,'ascend');

% show each events
tax_glu_sub=[-GluLag:GluLag]*GluTimeGap;
gspMat2=gspMat;
gspMat2(gspMat2==0)=NaN;

figure(11); clf; ax1=[];
for ev=1:size(gspMat,2)
    gspMat_sub=squeeze(gspMat2(:,ev,:));
    gspMat_sub_sorted=squeeze(gspMat2(GluROIsorted,ev,:));
    [rr, cc] = find(~isnan(gspMat_sub_sorted));
    ax1=[ax1 nexttile([1 1])];
    scatter(tax_glu_sub(cc)', rr, 20, 'k', '|', 'LineWidth', 2); hold all;
    scatter(tax_glu_sub,gspMat_sub(Syn2show,:),'r|','LineWidth',3); 
    ylim([-1 50])
    yyaxis right;
    plot(tax_glu_sub,sum(~isnan(gspMat_sub),1),'b-');
    ylim([0 15])
    ax1=[ax1 nexttile([1 1])];
    Volt2show=squeeze(gspMat_voltage(voltROIsorted,ev,:));
    imagesc([-VoltLag:VoltLag]/1000,[1:nROI],movmean(Volt2show,15,2),caxis); colormap(turbo);
end
linkaxes(ax1,'x');
xlim([-150 150]/1000)
%%
% show voltage decay of cluster and isolated events
figure(12); clf; tiledlayout(2,4); 
pixelsize=0.4680; %µm
coincidence_threshold=1;
distancethreshold=50;
coactiveSyn_low_thresh=2; coactiveSyn_upper_thresh=5;

gspMat_coincidence=sum(gspMat(:,:,round(end/2)+[-coincidence_threshold:coincidence_threshold]),[3],'omitnan');
NcoactiveSyn=sum(gspMat_coincidence(GluFtprntDistFromSyn*pixelsize<distancethreshold,:));
EventIsolated=NcoactiveSyn<=coactiveSyn_low_thresh;
EventClustered=NcoactiveSyn>=coactiveSyn_upper_thresh;
IsolatedGTA=squeeze(mean(gspMat_voltage(voltROIsorted,EventIsolated,:),2));
%IsolatedGTA=IsolatedGTA-median(IsolatedGTA,2);
ClusteredGTA=squeeze(mean(gspMat_voltage(voltROIsorted,EventClustered,:),2));
%ClusteredGTA=ClusteredGTA-median(ClusteredGTA,2);

nexttile([1 1]);
histogram(NcoactiveSyn);
xlabel('# of coactive synapses'); ylabel('# of events'); box off;

nexttile([1 1]);
imagesc([-VoltLag:VoltLag]/1000,[1:nROI],IsolatedGTA,caxis); hold all; yyaxis right;
clusterActivity=squeeze(mean(gspMat(find(GluFtprntDistFromSyn*pixelsize<distancethreshold),find(EventIsolated),:),[1 2],'omitnan'));
plot(tax_glu_sub,clusterActivity,'linewidth',2); ylabel(sprintf('Fraction of coactive synapses within %d µm',distancethreshold));
title(sprintf(['Isolated glutamate-event \n triggered average voltage \n,' ...
    ' Coactive synapse ≤ %d within %d µm)'],coactiveSyn_low_thresh,distancethreshold));
xlim([-150 150]/1000); ylim([0 0.15]);
xlabel('Time (sec)');

nexttile([1 1]);
imagesc([-VoltLag:VoltLag]/1000,[1:nROI],ClusteredGTA,caxis); hold all; yyaxis right;
title(sprintf(['Clustered glutamate-event \n triggered average voltage \n,' ...
    ' Coactive synapse ≥ %d within %d µm)'],coactiveSyn_upper_thresh,distancethreshold));
clusterActivity=squeeze(mean(gspMat(find(GluFtprntDistFromSyn*pixelsize<distancethreshold),find(EventClustered),:),[1 2],'omitnan'));
plot(tax_glu_sub,clusterActivity,'linewidth',2); ylim([0 0.15]);
ylabel(sprintf('Fraction of coactive synapses within %d µm',distancethreshold))
xlabel('Time (sec)');
xlim([-150 150]/1000); 

nexttile([1 1]);
x_fit=[min(InterGludist_sorted):max(InterGludist_sorted)]'*pixelsize;
[y_fit_iso, t_consts_iso] = expfitDM_2(InterGludist_sorted'*pixelsize, mean(IsolatedGTA(:,VoltLag+[-30:30]),2),x_fit, 100);
[y_fit_clu, t_consts_clu] = expfitDM_2(InterGludist_sorted'*pixelsize, mean(ClusteredGTA(:,VoltLag+[-30:30]),2),x_fit, 100);

plot(InterGludist_sorted*pixelsize,mean(IsolatedGTA(:,VoltLag+[-30:30]),2),'color',[0 0.2 1]); hold all
plot(InterGludist_sorted*pixelsize,mean(ClusteredGTA(:,VoltLag+[-30:30]),2),'color',[1 0.2 0]);
plot(x_fit,y_fit_iso,'b--');
plot(x_fit,y_fit_clu,'r--');
legend({sprintf('Isolated, Coactive synapse \n within %d µm ≤ %d, n = %d event',distancethreshold,coactiveSyn_low_thresh,sum(EventIsolated)),...
        sprintf('clustered, Coactive synapse \n within %d µm ≥ %d, n = %d event',distancethreshold,coactiveSyn_upper_thresh,sum(EventClustered)),...
        sprintf('Isolated fit, \\it{L} = %3.0f µm',t_consts_iso),...
        sprintf('Clusterd fit, \\it{L} = %3.0f µm',t_consts_clu)})
xlabel('Contour distance from synapse of interest (µm)');
ylabel('Voltage (-∆F/F)'); box off;

nexttile([1 4]);
show_footprnt_contour(AlignedVolt_ftprnt(:,:,voltROIsorted),VoltRefIm); hold all
scatter(Glu_coord(Syn2show,1),Glu_coord(Syn2show,2),'ro')
drawScaleBar(100/pixelsize,'horizontal')

%% Cluster size and subthreshold
tau=1; dist_threshold=30; %µm
GluMatched2Neuron= matchSglu2SWC(dendritePoints, Glu_coord, 'MaxDist', 5);
Syn_N1=GluMatched2Neuron.neuronID==1;
GluSpike_sub=GluSpike(Syn_N1,:);
Glu_coord_sub=Glu_coord(Syn_N1,:);
nSyn=sum(Syn_N1);
SkelDend = DendriteMask_N1;

%load SynDistance_sub
load('/Volumes/BHL18TB_D2/20260203_SD_V2+iGluSNFR4/170544BHLm302_V_Glu_25X_3min_N2/DistanceMatrix_GluSignal_N1.mat')
SynDistance_sub=[];
for i=1:nSyn
    for j=i:nSyn
    [SynDistance_sub(i,j), p]=geodesic_distance(SkelDend,Glu_coord_sub(i,:),Glu_coord_sub(j,:));
    end
    if mod(i, 20) == 0
        fprintf('  Frame %d / %d  (%.1f%%)\n', i, nSyn, 100*i/nSyn);
    end
end
load('/Volumes/BHL18TB_D2/20260203_SD_V2+iGluSNFR4/170544BHLm302_V_Glu_25X_3min_N2/DistanceMatrix_GluSignal_N1.mat')

GluEvents_cluster= find_spatiotemporal_events(GluSpike_sub, SynDistance_sub, Glu_coord_sub,'MaxDist',10/pixelsize,'Tolerance',1);
EventCentroid=vertcat(GluEvents_cluster.centroid);
CentroidFtprnt_distMat=distance_mat(EventCentroid,get_coord(AlignedVolt_ftprnt));
[Cent2Ftprnt_dist, minarg]=min(CentroidFtprnt_distMat,[],2);
validEvent=Cent2Ftprnt_dist*pixelsize<dist_threshold;
GluEvents_cluster2=GluEvents_cluster(find(validEvent));
minarg=minarg(validEvent);
GluFtprnt_distMat=distance_mat(Glu_coord(Syn_N1,:),get_coord(AlignedVolt_ftprnt));
[Glu2Ftprnt_dist, minarg2]=min(GluFtprnt_distMat,[],2);

[eventinVolt] = match_nearest(GluResult.t_ax(vertcat(GluEvents_cluster(find(validEvent)).time)), VoltResult.tax);
[~, eventMat_volt, eventmatchedVoltageTime]=get_STA(VoltTrace_Norm,eventinVolt,VoltLag,VoltLag);
eventMat_volt=eventMat_volt-median(eventMat_volt,3);
isVoltValid=ismember(eventinVolt,eventmatchedVoltageTime);
GluEvents_voltValid=GluEvents_cluster2(isVoltValid);

EventclusterSize=vertcat(GluEvents_voltValid.size);
minarg=minarg(isVoltValid);
EventVoltage=[];
for ev=1:size(eventMat_volt,2)
EventVoltage(ev,:)=eventMat_volt(minarg(ev),ev,:);
end
EventVoltageClustered=[];
for clsz=[min(EventclusterSize):max(EventclusterSize)]
    EventVoltageClustered{clsz}=mean(EventVoltage(EventclusterSize==clsz,VoltLag+[-30:30]),2,'omitnan');
end

figure(13); Violin_wPoints(EventVoltageClustered,hsv(9))
set(gca,'XTickLabel',[{'Isolated'}, arrayfun(@num2str, 2:9, 'UniformOutput', false)]);
xlabel('Cluster size (# synapse)'); ylabel('Subthreshold at glutamate event (∆F/F)');

% Depolarization variation across synapses at isolated event 
roi_event_idx  = cell(nSyn, 1);    % event indices for each ROI
roi_event_size = cell(nSyn, 1);    % event sizes   for each ROI
for n = 1:nSyn
    roi_event_idx{n}  = find(arrayfun(@(e) ismember(n, e.rois), GluEvents_voltValid));
    roi_event_size{n} = [GluEvents_voltValid(roi_event_idx{n}).size];
end

clustersize2show=[1];
median_clustersize=cellfun(@(x) prctile(x,90),roi_event_size);
idx2show=cellfun(@(x,y) y(ismember(x,clustersize2show)),roi_event_size,roi_event_idx,'UniformOutput', false);
IsolatedSubthreshold=[];
for n=1:nSyn   
subthvec=mean(mean(eventMat_volt(minarg2(n),idx2show{n},VoltLag+[-30:30]),[3]),2);
IsolatedSubthreshold=[IsolatedSubthreshold; [n median_clustersize(n) subthvec]];
end

figure(14); clf; tiledlayout(1,2);
ninetyprc_threshold=5;
nexttile([1 1]);
scatter(IsolatedSubthreshold(:,2),IsolatedSubthreshold(:,3),10,'k','filled')
xlabel('90th percentile of cluster size'); ylabel(sprintf('Subthreshold at isolated \n glutamate event (∆F/F)'));
nexttile([1 1]);
IsolatedlikeSyn=IsolatedSubthreshold(:,2)<ninetyprc_threshold;
ClusterlikeSyn=IsolatedSubthreshold(:,2)>ninetyprc_threshold;
p=Boxplot_wPoints2({IsolatedSubthreshold(IsolatedlikeSyn,3) IsolatedSubthreshold(ClusterlikeSyn,3)},hsv(2));
drawPValueLines(p,0,'StepHeight',0.005,'TextYOffset',0.002); box off;
ylim([-8 14]*0.001); ylabel(sprintf('Subthreshold at isolated \n glutamate event (∆F/F)'))
set(gca,'XTickLabel',{'Isolated-like synapse','Cluster-like synapse'})

%% Glu-triggered voltage movie
% For each glutamate event time, extract N frames before and M frames after
% from the voltage movie, apply all artefact regressions, average across
% events, and save as a 3-panel video (Glu frame | single-event Volt | STA Volt).
%
% REQUIRES in workspace (from make_GluVolt_movie setup block):
%   fpath, f, Device_Data, VoltResult, GluResult, tformReg
%   framesPerSeg, Fs_volt, volt_smooth_sigma
%   mov_glu_dFF_masked, t_glu_show   (processed Glu movie + its time axis)
%   cmap_volt, cmap_glu, volt_dFF_range, glu_dFF_range
%   DendriteMask, GluResult.AvgGluImg

% ===================== USER SETTINGS ====================================

% Trigger times: glutamate event times in SECONDS (same scale as t_glu_show)
% Example: use peak frames from GluResult or manual selection
Syn2look=7;
GluSpikeTime2look=find(GluSpike(Syn2look,:)>0);
glu_event_times_s = GluResult.t_ax(find(GluSpike(Syn2look,:)>0));

N_pre  = 500;    % voltage frames BEFORE each trigger  (ms at 1000 Hz)
M_post = 500;    % voltage frames AFTER  each trigger  (ms at 1000 Hz)

output_STA_filename = fullfile(fpath{f}, 'GluTriggered_VoltSTA.mp4');
frame_rate          = 30;
quality             = 95;

pixelsize=0.4680; %µm
coincidence_threshold=1;
distancethreshold=50;
coactiveSyn_low_thresh=3; coactiveSyn_upper_thresh=6;

nWin      = N_pre + M_post +1;          % total window length in volt frames
tau_axis  = (-N_pre : M_post);     % frame offset axis (0 = trigger)
tau_axis_glu= -ceil(N_pre/ (exposuretime2*1000)):ceil(M_post/ (exposuretime2*1000));

[glu_event_volt_frames , dist] = match_nearest(glu_event_times_s , VoltResult.t_ax);
[~, gspMat]=get_STA(GluSpike,GluSpike(Syn2show,:)>0,GluLag,GluLag);
[~, gspMat_voltage, validglu]=get_STA(VoltTrace_Norm,glu_event_volt_frames,N_pre,M_post);
validglu=ismember(gsptime_volt,validglu);
GluSpikeTime2look=GluSpikeTime2look(validglu);
glu_event_volt_frames=glu_event_volt_frames(validglu);
gspMat=gspMat(:,validglu,:);
gspMat_voltage=gspMat_voltage-median(gspMat_voltage,3);

gspMat_coincidence=sum(gspMat(:,:,round(end/2)+[-coincidence_threshold:coincidence_threshold]),[3],'omitnan');
NcoactiveSyn=sum(gspMat_coincidence(GluFtprntDistFromSyn*pixelsize<distancethreshold,:));
EventIsolated=NcoactiveSyn<=coactiveSyn_low_thresh;
EventClustered=NcoactiveSyn>=coactiveSyn_upper_thresh;
IsolatedGTA=squeeze(mean(gspMat_voltage(voltROIsorted,EventIsolated,:),2));
ClusteredGTA=squeeze(mean(gspMat_voltage(voltROIsorted,EventClustered,:),2));

%% ===================== PRE-COMPUTE VOLT REGRESSION REGRESSORS ===========
% These are the same regressors used in make_GluVolt_movie voltage block.
% We build them globally so we can slice the right segment per window.
fprintf('Preparing voltage event movies...\n');

% Exponential bleaching fit — compute on a reference segment (seg 1)
volt_smooth_sigma=1;
t2bleachfit=[15001:25001];
mov_ref_seg = readBinMov_BHL_multiple(fpath{f}, 3, t2bleachfit, framesPerSeg, 'mc');
meanF_ref   = squeeze(mean(mov_ref_seg, [1 2]));
y_fit_ref   = expfitDM_2(t2bleachfit', meanF_ref, (1:length(VoltResult.t_ax)-1)', 10000);
% Normalise so it can be scaled to any segment
bleach_ref  = y_fit_ref / mean(y_fit_ref);   % unit-mean bleach profile
mc_win = VoltResult.mc(t2bleachfit, :);
mov_ref_seg = SeeResiduals(mov_ref_seg, bleach_ref(t2bleachfit));
mov_ref_seg = SeeResiduals(mov_ref_seg, mc_win);
mov_ref_seg = SeeResiduals(mov_ref_seg, mc_win.^2);
mov_ref_seg = SeeResiduals(mov_ref_seg, mc_win(:,1) .* mc_win(:,end));
mov_ref_seg = SeeResiduals(mov_ref_seg, VoltResult.bvTrace(:, t2bleachfit), 1);

Fstd_img = get_F0img_PCA(imgaussfilt(mov_ref_seg(:,:,[1:5:10000]),volt_smooth_sigma));
F0_img = imgaussfilt(VoltResult.ref_im, volt_smooth_sigma);
clear mov_ref_seg

sz_volt  = size(VoltResult.ref_im);
sz_glu = size(GluResult.AvgGluImg);
nEvents= length(glu_event_volt_frames);

maskIdx=find(DendriteMask_N1>0);

TA_volt     = zeros(length(maskIdx), nWin, nEvents, 'double');
TA_volt_sub = zeros(length(maskIdx), nWin, nEvents, 'double');

for ei = 1:nEvents

    win_fr  = glu_event_volt_frames(ei) + tau_axis;          % global volt frames in window
    mov_ev = readBinMov_BHL_multiple(fpath{f}, 3, win_fr, framesPerSeg, 'mc');

    mov_ev_win = mov_ev - mean(mov_ev, 3);
    bkg_win     = zeros(1, nWin);
    bkg_win(1,:) = bleach_ref(win_fr);
    bkg_win = [bkg_win; VoltResult.bvTrace(:, win_fr)];

    mc_win = VoltResult.mc(win_fr, :);
    mov_ev_win = SeeResiduals(mov_ev_win, mc_win);
    mov_ev_win = SeeResiduals(mov_ev_win, mc_win.^2);
    mov_ev_win = SeeResiduals(mov_ev_win, mc_win(:,1) .* mc_win(:,end));
    mov_ev_win = SeeResiduals(mov_ev_win, bkg_win, 1);
    
    dFF_win     = -tovec(imgaussfilt(mov_ev_win, volt_smooth_sigma) ./ Fstd_img);
    %dFF_win     = -tovec(imgaussfilt(mov_ev_win, volt_smooth_sigma) ./ F0_img);
    dFF_win_sub = movmean(dFF_win, 20, 2);

    TA_volt(:,:,ei)=dFF_win(maskIdx,:);
    TA_volt_sub(:,:,ei)=dFF_win_sub(maskIdx,:);
    fprintf('  Event %d / %d done.\n', ei, nEvents);
end

%-- Normalise STA
STA_volt = zeros(sz_volt(1)*sz_volt(2),nWin);
STA_volt_sub = zeros(sz_volt(1)*sz_volt(2),nWin);
STA_volt_sub_Isolated  = zeros(sz_volt(1)*sz_volt(2),nWin);
STA_volt_sub_Clustered = zeros(sz_volt(1)*sz_volt(2),nWin);

STA_volt(maskIdx,:)     = mean(TA_volt,3,'omitnan');
STA_volt_sub(maskIdx,:) = mean(TA_volt_sub,3,'omitnan');
STA_volt_sub_Isolated(maskIdx,:)  = mean(TA_volt_sub(:,:,EventIsolated),3,'omitnan');
STA_volt_sub_Clustered(maskIdx,:) = mean(TA_volt_sub(:,:,EventClustered),3,'omitnan');

STA_volt=toimg(STA_volt,sz_volt(1),sz_volt(2));
STA_volt_sub=toimg(STA_volt_sub,sz_volt(1),sz_volt(2));
STA_volt_sub_Isolated=toimg(STA_volt_sub_Isolated,sz_volt(1),sz_volt(2));
STA_volt_sub_Clustered=toimg(STA_volt_sub_Clustered,sz_volt(1),sz_volt(2));

[~, TA_Glu]=get_STA(tovec(dFF_roi),GluSpikeTime2look,-tau_axis_glu(1),tau_axis_glu(end));
STA_Glu=toimg(squeeze(mean(TA_Glu(:,:,:),2,'omitnan')),sz_glu(1),sz_glu(2));
STA_Glu_Isolated=toimg(squeeze(mean(TA_Glu(:,EventIsolated,:),2,'omitnan')),sz_glu(1),sz_glu(2));
STA_Glu_clustered=toimg(squeeze(mean(TA_Glu(:,EventClustered,:),2,'omitnan')),sz_glu(1),sz_glu(2));

TAxis_GluTrigger=tau_axis_glu*exposuretime2;
TAxis_VoltTrigger=tau_axis*exposuretime1;
STA_Glu_interp = interp1(TAxis_GluTrigger, tovec(STA_Glu)', TAxis_VoltTrigger, 'linear', 'extrap').';
STA_Glu_Isolated_interp = interp1(TAxis_GluTrigger, tovec(STA_Glu_Isolated)', TAxis_VoltTrigger, 'linear', 'extrap').';
STA_Glu_clustered_interp = interp1(TAxis_GluTrigger, tovec(STA_Glu_clustered)', TAxis_VoltTrigger, 'linear', 'extrap').';

STA_Glu_interp =toimg(STA_Glu_interp ,sz_glu(1),sz_glu(2));
STA_Glu_Isolated_interp =toimg(STA_Glu_Isolated_interp ,sz_glu(1),sz_glu(2));
STA_Glu_clustered_interp =toimg(STA_Glu_clustered_interp ,sz_glu(1),sz_glu(2));

STA_Glu_interp=STA_Glu_interp-median(STA_Glu_interp,3);
STA_Glu_Isolated_interp=STA_Glu_Isolated_interp-median(STA_Glu_Isolated_interp(:,:,1:200),3);
STA_Glu_clustered_interp=STA_Glu_clustered_interp-median(STA_Glu_clustered_interp(:,:,1:200),3);

STA_Glu_interp(isnan(STA_Glu_interp))=0;
STA_Glu_Isolated_interp(isnan(STA_Glu_Isolated_interp))=0;
STA_Glu_clustered_interp(isnan(STA_Glu_clustered_interp))=0;

%% Polyline kymo of STA movie

Vmeta = Device_Data{1,3};
Gmeta = Device_Data{1,4};
SynCoord2show=Glu_coord(Syn2show,:);
SynCoord2show_volt=zeros(size(SynCoord2show));
Rvolt = refFromROI(size(VoltResult.ref_im), double(Vmeta.ROI));
Rglu  = refFromROI(size(GluResult.AvgGluImg),double(Gmeta.ROI));
[xWorld, yWorld] = intrinsicToWorld(Rglu, SynCoord2show(:,1), SynCoord2show(:,2));
[SynCoord2show_volt(:,1), SynCoord2show_volt(:,2)] = transformPointsForward(tformReg, xWorld, yWorld);
[SynCoord2show_volt(:,1), SynCoord2show_volt(:,2)] = worldToIntrinsic(Rvolt, SynCoord2show_volt(:,1), SynCoord2show_volt(:,2));

[KymoTr KymoROI]=polyLineKymo3(STA_volt_sub_Clustered,10,10,VoltResult.ref_im);
KymoROI_coord=cellfun(@(x) mean(x(1:end-1,:)),KymoROI,'UniformOutput',false);
KymoROI_coord=cell2mat(KymoROI_coord');

Depol_glu=mean(KymoTr(500+[-15:5],:),1);
Depol_cmap=vec2cmap(Depol_glu,turbo(256));
ContourDist=[];
for roi=1:length(KymoROI)
ContourDist(roi)=geodesic_distance(DendriteMask_N1,KymoROI_coord(roi,:),SynCoord2show_volt);
end
[~, minROI]=min(ContourDist);

figure(22); clf;
tiledlayout(2,2);
nexttile([1 2]);
imshow2(VoltResult.ref_im,[]); hold all;
scatter(KymoROI_coord(:,1),KymoROI_coord(:,2),25,Depol_cmap,'filled');
l=scatter(SynCoord2show_volt(:,1),SynCoord2show_volt(:,2),50,'g^','filled'); hold all;
legend(l,'Triggering synapse')

nexttile([1 1]);
imagesc([tau_axis],[1:length(KymoROI)],KymoTr'); hold all;
scatter([tau_axis(round(end/2))-145],minROI,40,'r>','filled');
xlim([-150 150]); xlabel('Time (ms)'); ylabel('ROI ID');

nexttile([1 1]);
plot(ContourDist(1:minROI)*pixelsize,Depol_glu(1:minROI)); hold all
plot(ContourDist(minROI:end)*pixelsize,Depol_glu(minROI:end));
xlabel('Contour distance (µm)'); ylabel('Subthreshold (∆F/F)'); box off;
legend({'Tip to synapse','Synapse to soma'})



%% ===================== VIDEO WRITER =====================================
v = VideoWriter(output_STA_filename, 'MPEG-4');
v.FrameRate = frame_rate;
v.Quality   = quality;
open(v);
fig = figure('Color','k','Units','normalized','Position',[0.2 0.2 0.6 0.5]);
fprintf('Rendering STA movie...\n');
vcaxis=[-1 1]; gcaxis=[20 100];
OverlayScale=[1 0.5];
Str_rgb=grs2rgb(VoltResult.ref_im,gray(256),20, 200);
bd=5;
cmap_volt = gen_colormap([0 0.2 1; 0 0 0; 1 0 0], 256);
cmap_glu  = gen_colormap([0 0 0; 0 1 0], 256);

Alpha_thres=0.5;
Alpha_steepness=10;
Alpha_kernal=4;

for fi = 1:2:length(tau_axis)

    alpha = 1 ./ (1 + exp(-Alpha_steepness * (abs(imgaussfilt(STA_volt_sub(bd:end-bd,bd:end-bd,fi),Alpha_kernal)) - Alpha_thres)));
    sub_rgb=grs2rgb(STA_volt_sub(bd:end-bd,bd:end-bd,fi),cmap_volt,vcaxis(1),vcaxis(2)).*alpha.*DendriteMask_N1(bd:end-bd,bd:end-bd);
    alpha2 = 1 ./ (1 + exp(-Alpha_steepness * (abs(imgaussfilt(STA_volt_sub_Isolated(bd:end-bd,bd:end-bd,fi),Alpha_kernal)) - Alpha_thres)));
    sub_iso_rgb=grs2rgb(STA_volt_sub_Isolated(bd:end-bd,bd:end-bd,fi),cmap_volt,vcaxis(1),vcaxis(2)).*alpha2.*DendriteMask_N1(bd:end-bd,bd:end-bd);
    alpha3 = 1 ./ (1 + exp(-Alpha_steepness * (abs(imgaussfilt(STA_volt_sub_Clustered(bd:end-bd,bd:end-bd,fi),Alpha_kernal)) - Alpha_thres)));
    sub_clu_rgb=grs2rgb(STA_volt_sub_Clustered(bd:end-bd,bd:end-bd,fi),cmap_volt,vcaxis(1),vcaxis(2)).*alpha3.*DendriteMask_N1(bd:end-bd,bd:end-bd);

    SGA_rgb= transformCamera(Device_Data, tformReg, VoltResult.ref_im, STA_Glu_interp(:,:,fi));
    SGA_iso_rgb= transformCamera(Device_Data, tformReg, VoltResult.ref_im, STA_Glu_Isolated_interp(:,:,fi));
    SGA_clu_rgb= transformCamera(Device_Data, tformReg, VoltResult.ref_im, STA_Glu_clustered_interp(:,:,fi));

    SGA_rgb=grs2rgb(SGA_rgb,cmap_glu,gcaxis(1),gcaxis(2));
    SGA_iso_rgb=grs2rgb(SGA_iso_rgb,cmap_glu,gcaxis(1),gcaxis(2));
    SGA_clu_rgb=grs2rgb(SGA_clu_rgb,cmap_glu,gcaxis(1),gcaxis(2));

    %-- Render
    clf(fig);
    set(fig, 'Color','k');
    tl = tiledlayout(fig, 3, 3, 'TileSpacing','compact', 'Padding','compact');

    ax7 = nexttile(tl, 1);
    imagesc(imgaussfilt(STA_volt_sub(bd:end-bd,bd:end-bd,fi),1),[vcaxis]); axis equal tight off;
    title(ax7, sprintf('Subthreshold voltage (%+.0f ms)', TAxis_VoltTrigger(fi)*1000), 'Color','w', 'FontSize',11);
    colormap(gray);

    ax7 = nexttile(tl, 4);
    imagesc(imgaussfilt(STA_volt_sub_Isolated(bd:end-bd,bd:end-bd,fi),1),[vcaxis]); axis equal tight off;
    colormap(gray);

    ax7 = nexttile(tl, 7);
    imagesc(imgaussfilt(STA_volt_sub_Clustered(bd:end-bd,bd:end-bd,fi),1),[vcaxis]); axis equal tight off;
    colormap(gray);

    ax1 = nexttile(tl, 2);
    imagesc(sub_rgb*OverlayScale(1)+Str_rgb(bd:end-bd,bd:end-bd,:)*OverlayScale(2)); axis equal tight off;
    title(ax1, sprintf('Average of all glutamate event'), 'Color','w', 'FontSize',11);

    ax2 = nexttile(tl, 3);
    imagesc(SGA_rgb(bd:end-bd,bd:end-bd,:)+Str_rgb(bd:end-bd,bd:end-bd,:)*OverlayScale(2)); axis equal tight off;
    title(ax2, sprintf('Glutamate'), 'Color','w', 'FontSize',11);

    ax3 = nexttile(tl, 5);
    imagesc(sub_iso_rgb*OverlayScale(1)+Str_rgb(bd:end-bd,bd:end-bd,:)*OverlayScale(2)); axis equal tight off;
    title(ax3, sprintf('Average of isolated glutamate event,\n N = %2.0d',sum(EventIsolated)), 'Color','w', 'FontSize',11);

    ax4 = nexttile(tl, 6);
    imagesc(SGA_iso_rgb(bd:end-bd,bd:end-bd,:)+Str_rgb(bd:end-bd,bd:end-bd,:)*OverlayScale(2)); axis equal tight off;

    ax5 = nexttile(tl, 8);
    imagesc(sub_clu_rgb*OverlayScale(1)+Str_rgb(bd:end-bd,bd:end-bd,:)*OverlayScale(2)); axis equal tight off;
    title(ax5, sprintf('Average of clustered glutamate event,\n N = %2.0d',sum(EventClustered)'), 'Color','w', 'FontSize',11);

    ax6 = nexttile(tl, 9);
    imagesc(SGA_clu_rgb(bd:end-bd,bd:end-bd,:)+Str_rgb(bd:end-bd,bd:end-bd,:)*OverlayScale(2)); axis equal tight off;
    drawScaleBar(100/pixelsize,'horizontal');

    drawnow;
    writeVideo(v, getframe(fig));

    if mod(fi, 50) == 0
        fprintf('  Frame %d / %d  (%.1f%%)\n', fi, nWin, 100*fi/nWin);
    end
end

close(v);
close(fig);
fprintf('Done. STA movie saved to:\n  %s\n', output_STA_filename);
