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

%% ===================== USER SETTINGS ====================================
f            = 3;        % session index into fpath{}

% Voltage frame range to render (global 1-based frame indices).
% Frame 1 = first frame of mc01.bin, frame 15001 = first frame of mc02.bin, etc.
t_start_volt = 12000;        % e.g.     1  ->  0.000 s  (at 1000 Hz)
t_end_volt   = 30000;     % e.g.  3000  ->  3.000 s

framesPerSeg  = 15000;   % frames per mc0x.bin file  (15 s × 1000 Hz)
Fs_volt       = 1000;    % voltage camera frame rate (Hz)

% Spatial smoothing (Gaussian sigma in pixels, per frame)
volt_smooth_sigma = 1.0;
glu_smooth_sigma  = 1.5;

% Output video
output_filename = fullfile(fpath{f}, 'GluVolt_movie.mp4');
frame_rate      = 30;      % fps of output video
quality         = 95;

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

%% ===================== BUILD VOLTAGE TIME AXIS ==========================
% Simple: frame index IS the time (no DAQ clock needed).
% Global frame 1 corresponds to t = 0 s.
t_start_volt = max(t_start_volt, 1);
t_show       = ([t_start_volt:t_end_volt]- 1) / Fs_volt; % seconds
frm2read_v = [t_start_volt:t_end_volt];
frm2read_g = find(GluResult.t_ax>=t_show(1) & GluResult.t_ax<t_show(end));

fprintf('Rendering %d voltage frames  (%.3f – %.3f s)\n', ...
    numel(t_show), t_show(1), t_show(end));

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
t_glu = GluResult.t_ax;   % cam2_trig is in volt-frame units

% Rebuild GluemovTimesegments for loading
if cam2_vsyn(end) < GlumoviemaxTime
    GluemovTimesegments = [cam2_vsyn(mod488)+2,  cam2_vsyn(end)];
else
    GluemovTimesegments = [(cam2_vsyn(mod488)+2) : GlumoviemaxTime : cam2_vsyn(end)];
    GluemovTimesegments(end+1) = cam2_vsyn(end);
end

%% ===================== PRELOAD GLU MOVIE ================================
fprintf('Loading glutamate movie segments...\n');
mov_glu =readBinMov_BHL_multiple(fpath{f},4,(frm2read_g),5000,'mc2');
szGlu= size(mov_glu);
[V, D, u_sub] = get_eigvector(tovec(imgaussfilt(mov_glu,glu_smooth_sigma))',40);
u_sub=reshape(u_sub,szGlu(1),szGlu(2),[]);
figure(2); clf; imshow2_patch(u_sub); drawnow;
n=input("PCs to regress out\n");
mov_glu2 =SeeResiduals(mov_glu ,V(:,n));
mov_glu2 =SeeResiduals(mov_glu2 ,GluResult.mc([frm2read_g],:));

% % Glu dF/F
% mov_glu2=imgaussfilt(mov_glu2,glu_smooth_sigma);
% glu_F0      = movmedian(mov_glu, 200, 3);
% mov_glu_dFF = (mov_glu2 - glu_F0) ./ (glu_F0 + 1);
% mov_glu_dFF = mov_glu_dFF-median(mov_glu_dFF,3);
% glu_ref_im  = mean(mov_glu, 3);
% GluMask=(max(GluResult.S_glu>0,[],3)>0);
% GluMask= imdilate(GluMask, strel('disk', 2));
% mov_glu_dFF_masked=mov_glu_dFF.*(max(GluResult.S_glu>0,[],3)>0);
% 
% % Reshape to (T x XY), interpolate, reshape back
% mov2d       = tovec(mov_glu_dFF_masked)';          % T  x XY
% mov2d_interp = interp1(t_glu(frm2read_g), mov2d, t_show, 'linear','extrap');  % T2 x XY
% mov_glu_dFF_interp= reshape(mov2d_interp', size(mov_glu_dFF,1), size(mov_glu_dFF,2), []); % X x Y x T2

% Glu dF/F — per-ROI masked, median-filtered, and normalized
n_kernel = 2.5;   % median filter kernel size in pixels (set as needed)
dil_rad  = 2;   % dilation radius in pixels

N_roi = size(GluResult.S_glu, 3);

% Pre-allocate
mov_glu_dFF_masked = zeros(size(mov_glu2), 'double');

%-- Normalize by median of raw mov_glu within this ROI
    % Compute per-pixel mean over time within the (undilated) ROI
    F0_roi=tovec(mean(mov_glu,3))'* tovec(GluResult.S_glu)./squeeze(sum(tovec(GluResult.S_glu),[1 2]));

for i = 1:N_roi

    %-- (i) Build dilated mask for this ROI
    roi_mask     = GluResult.S_glu(:,:,i) > 0;               % X x Y binary
    roi_mask_dil = imdilate(roi_mask, strel('disk', dil_rad));% dilated mask

    %-- Extract movie pixels within dilated mask, set outside to NaN
    % Apply Gaussian smooth first (same as before), then mask
    mov_roi = mov_glu2;               
    mov_roi=maskBinary(mov_roi,roi_mask_dil==0,NaN);

    %-- (ii) 2D median filter within masked area, per frame
    mov_roi= imgaussfiltnan(mov_roi, n_kernel); 

    % dF/F per pixel: (filtered - F0) / F0
    dFF_roi = (mov_roi) ./ (F0_roi(i));

    % Zero-centre per pixel by subtracting median over time
    %dFF_roi = dFF_roi - median(dFF_roi, 3, 'omitnan');

    % Accumulate into full movie (only within dilated mask)
    for t = 1:size(dFF_roi, 3)
        frame = dFF_roi(:,:,t);
        frame(~roi_mask_dil) = 0;
        mov_glu_dFF_masked(:,:,t) = mov_glu_dFF_masked(:,:,t) + frame;
    end

    if mod(i, 10) == 0
        fprintf('\r  ROI %d / %d  (%.0f%% complete)', i, N_roi, 100*i/N_roi);
        drawnow;
    end
end
fprintf('\n');

%-- Glu time axis for per-frame interpolation in render loop
t_glu_show = t_glu(frm2read_g);   % time axis of mov_glu_dFF_masked (T_glu x 1)

%% ===================== VOLTAGE F0 IMAGE =================================
% Load first complete segment to build a pixel-wise F0
fprintf('Building voltage F0 image from mc01.bin...\n');
mov_volt=readBinMov_BHL_multiple(fpath{f},3,frm2read_v,framesPerSeg,'mc');
mov_res= mov_volt-mean(mov_volt,3);
meanF=squeeze(mean(mov_volt,[1 2]));
y_fit=expfitDM_2([1:size(mov_res,3)]',meanF,[1:size(mov_res,3)]',10000);
bkg = zeros(1, size(mov_volt,3));
bkg(1,:) =y_fit;
bkg =[bkg; VoltResult.bvTrace(:,frm2read_v)];
mov_res = SeeResiduals(mov_res,VoltResult.mc(frm2read_v,:));
mov_res = SeeResiduals(mov_res,VoltResult.mc(frm2read_v,:).^2);
mov_res = SeeResiduals(mov_res,VoltResult.mc(frm2read_v,1).*VoltResult.mc(frm2read_v,end));
mov_res= SeeResiduals(mov_res,bkg,1);
dFF_mov_volt=-imgaussfilt(mov_res,volt_smooth_sigma)./imgaussfilt(VoltResult.ref_im,volt_smooth_sigma);
dFF_mov_volt_sub=movmean(dFF_mov_volt,20,3);

%% ===================== COLORMAPS ========================================
volt_dFF_range = [-0.22  0.22];
glu_dFF_range  = [50  80];

cmap_volt = gen_colormap([0 0.2 1; 0 0 0; 1 0 0],256);

cmap_glu  = gen_colormap([0 0 0; 0 1 0], 256);

[H_out, W_out] = size(mov_glu);

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
DendriteMask=point2img(cell2mat(dendritePoints'),3,size(GluResult.AvgGluImg));
DendriteMask_N1=point2img(cell2mat(dendritePoints(1)),2,size(GluResult.AvgGluImg));

VoltRefIm= transformCamera_O2B(Device_Data, tformReg, VoltResult.ref_im, GluResult.AvgGluImg);
VoltRefIm_rgb = grs2rgb(VoltRefIm, gray(256), 20, 200);

%% ===================== VIDEO WRITER =====================================
% Colourmap display ranges

v = VideoWriter(output_filename, 'MPEG-4');
v.FrameRate = frame_rate;
v.Quality   = quality;
open(v);
fig = figure('Color','k', 'units','pixels', 'Position',[200 0 1000 800]);

% ===================== RENDER LOOP ======================================
fprintf('Rendering frames...\n');

for fi = 1:numel(frm2read_v)
    
    alignedVolt = transformCamera_O2B(Device_Data, tformReg, dFF_mov_volt(:,:,fi), GluResult.AvgGluImg);
    volt_rgb = grs2rgb(alignedVolt , cmap_volt, volt_dFF_range(1), volt_dFF_range(2));
    volt_rgb = volt_rgb .* DendriteMask .* (abs(alignedVolt)>0);

    alignedSub = transformCamera_O2B(Device_Data, tformReg, dFF_mov_volt_sub(:,:,fi), GluResult.AvgGluImg);
    alignedSub = maskBinary(alignedSub,DendriteMask==0,NaN);
    alignedSub = imgaussfiltnan(alignedSub ,2);

    voltSub_rgb = grs2rgb(alignedSub , cmap_volt, volt_dFF_range(1), volt_dFF_range(2));
    voltSub_rgb = voltSub_rgb .* DendriteMask .* (abs(alignedVolt)>0);

    %-- Per-frame Glu interpolation to current volt time
    gi_frac   = interp1(t_glu_show, 1:numel(t_glu_show), t_show(fi), 'linear', 'extrap');
    gi_frac   = min(max(gi_frac, 1), size(mov_glu_dFF_masked, 3));
    gi_lo     = max(floor(gi_frac), 1);
    gi_hi     = min(ceil(gi_frac),  size(mov_glu_dFF_masked, 3));
    alpha     = gi_frac - gi_lo;
    glu_frame = (1-alpha) * mov_glu_dFF_masked(:,:,gi_lo) + alpha * mov_glu_dFF_masked(:,:,gi_hi);
    glu_rgb   = grs2rgb(glu_frame, cmap_glu, glu_dFF_range(1), glu_dFF_range(2));

    %-- (G) Two-colour overlay: Glu=Green, Volt=Red
    overlay = zeros(H_out, W_out, 3);
    overlay = glu_rgb*1 + voltSub_rgb*1; %+ VoltRefIm_rgb*0.3;

    %-- (H) Render 3-panel figure
    clf(fig);
    set(fig, 'Color','k');
    tl = tiledlayout(fig, 2, 2, 'TileSpacing','tight', 'Padding','tight');

    ax1 = nexttile(tl, 1);
    imshow(glu_rgb, 'Parent', ax1);
    title(ax1, sprintf('Glutamate  t = %.3f s', t_show(fi)), 'Color','w', 'FontSize',11);

    ax2 = nexttile(tl, 2);
    imshow(volt_rgb, 'Parent', ax2);
    title(ax2, sprintf('Voltage  frame %d', frm2read_v(fi)), 'Color','w', 'FontSize',11);

    ax3 = nexttile(tl, 3);
    imshow(voltSub_rgb, 'Parent', ax3);
    title(ax2, sprintf('Voltage  frame %d', frm2read_v(fi)), 'Color','w', 'FontSize',11);

    ax4 = nexttile(tl, 4);
    imshow(overlay, 'Parent', ax4);
    title(ax3, 'Overlay  (G=Glu, R=Volt)', 'Color','w', 'FontSize',11);

    drawnow;
    writeVideo(v, getframe(fig));

    if mod(fi, 100) == 0
        fprintf('  Frame %d / %d  (%.1f%%)\n', fi, numel(frm2read_v), ...
            100*fi/numel(frm2read_v));
    end
end

close(v);
close(fig);
fprintf('Done. Movie saved to:\n  %s\n', output_filename);

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
[gsptime_volt, dist] = match_nearest(gsp_time, VoltResult.tax);
[gspTA_voltage(:,:,si), gspMat_voltage]=get_STA(VoltTrace_Norm,gsptime_volt,100,100);
end
gspTA_voltage=gspTA_voltage-median(gspTA_voltage,2);
Glu_coord=get_coord(GluResult.S_glu);

%% Show representative synaptic event case.

VoltLag=1000; %ms
Syn2show=309;
caxis=[-0.005 0.01];
GluTimeGap=min(diff(GluResult.t_ax));
GluLag=ceil(VoltLag/(GluTimeGap*1000));
gsp_time=GluResult.t_ax(find(GluSpike(Syn2show,:)>0));
[gsptime_volt, dist] = match_nearest(gsp_time, VoltResult.tax);
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

% show voltage decay of cluster and isolated events
figure(12); clf; tiledlayout(2,4); 
pixelsize=0.4680; %µm
coincidence_threshold=1;
distancethreshold=50;
coactiveSyn_low_thresh=3; coactiveSyn_upper_thresh=6;

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
% SynDistance_sub=[];
% for i=1:nSyn
%     for j=i:nSyn
%     [SynDistance_sub(i,j), p]=geodesic_distance(SkelDend,Glu_coord_sub(i,:),Glu_coord_sub(j,:));
%     end
%     if mod(i, 20) == 0
%         fprintf('  Frame %d / %d  (%.1f%%)\n', i, nSyn, 100*i/nSyn);
%     end
% end

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
glu_event_times_s = GluResult.t_ax(find(GluSpike(Syn2look,:)>0));

N_pre  = 250;    % voltage frames BEFORE each trigger  (ms at 1000 Hz)
M_post = 250;    % voltage frames AFTER  each trigger  (ms at 1000 Hz)

output_STA_filename = fullfile(fpath{f}, 'GluTriggered_VoltSTA.mp4');
frame_rate          = 30;
quality             = 95;

%% ===================== CONVERT GLU TIMES -> VOLT FRAME INDICES ==========
% t_glu_show is in seconds; volt frames are 1-based at Fs_volt = 1000 Hz

nWin      = N_pre + M_post;          % total window length in volt frames
tau_axis  = (-N_pre : M_post-1);     % frame offset axis (0 = trigger)
t_axis_ms = tau_axis;                % ms at 1000 Hz

[glu_event_volt_frames , dist] = match_nearest(gsp_time, VoltResult.tax);

% Remove events too close to edges
valid = glu_event_volt_frames > N_pre & ...
        glu_event_volt_frames <= (framesPerSeg * 99) - M_post;  % 99 = safe upper bound
glu_event_volt_frames = glu_event_volt_frames(valid);
nEvents = numel(glu_event_volt_frames);
fprintf('Using %d / %d events after edge exclusion.\n', nEvents, numel(glu_event_times_s));

%% ===================== PRE-COMPUTE VOLT REGRESSION REGRESSORS ===========
% These are the same regressors used in make_GluVolt_movie voltage block.
% We build them globally so we can slice the right segment per window.
fprintf('Preparing voltage regressors...\n');

% Exponential bleaching fit — compute on a reference segment (seg 1)
mov_ref_seg = readBinMov_BHL_multiple(fpath{f}, 3, 1:framesPerSeg, framesPerSeg, 'mc');
meanF_ref   = squeeze(mean(mov_ref_seg, [1 2]));
y_fit_ref   = expfitDM_2((1:framesPerSeg)', meanF_ref, (1:framesPerSeg)', 10000);
% Normalise so it can be scaled to any segment
bleach_ref  = y_fit_ref / mean(y_fit_ref);   % unit-mean bleach profile
clear mov_ref_seg

%% ===================== STA ACCUMULATOR ==================================
sz_volt  = [double(Device_Data{3}.ROI(4)), double(Device_Data{3}.ROI(2))];
nRow1    = sz_volt(1);  nCol1 = sz_volt(2);

STA_volt     = zeros(nRow1, nCol1, nWin, 'double');
STA_volt_sub = zeros(nRow1, nCol1, nWin, 'double');
nContrib     = zeros(1, nWin);

%% ===================== LOAD + PROCESS EACH EVENT ========================
cached_seg_idx  = -1;
cached_seg_data = [];

for ei = 1:nEvents

    trig_fr = glu_event_volt_frames(ei);   % global volt frame of trigger
    win_fr  = trig_fr + tau_axis;          % global volt frames in window

    %-- Identify which segments the window spans
    seg_indices = unique(ceil(win_fr / framesPerSeg));

    %-- Collect raw frames across (possibly multiple) segments
    mov_win = zeros(nRow1, nCol1, nWin, 'single');

    for fi = 1:nWin
        g  = win_fr(fi);
        si = ceil(g / framesPerSeg);           % bin file index
        lf = mod(g - 1, framesPerSeg) + 1;    % local frame within file

        if si ~= cached_seg_idx
            binfile = [fpath{f} '/mc' num2str(si,'%02d') '.bin'];
            cached_seg_data = single(readBinMov(binfile, nRow1, nCol1));
            cached_seg_idx  = si;
            fprintf('  Loaded %s\n', binfile);
        end
        mov_win(:,:,fi) = cached_seg_data(:,:,lf);
    end
    mov_win = double(mov_win);

    %-- Artefact regression (same pipeline as make_GluVolt_movie)
    % 1. Remove mean
    mov_res_win = mov_win - mean(mov_win, 3);

    % 2. Exponential bleaching (scaled bleach profile for this window)
    meanF_win   = squeeze(mean(mov_win, [1 2]));
    y_fit_win   = expfitDM_2((1:nWin)', meanF_win, (1:nWin)', 10000);
    bkg_win     = zeros(1, nWin);
    bkg_win(1,:) = y_fit_win;

    % 3. Blood vessel trace
    bkg_win = [bkg_win; VoltResult.bvTrace(:, win_fr)];

    % 4. Motion traces
    mc_win = VoltResult.mc(win_fr, :);
    mov_res_win = SeeResiduals(mov_res_win, mc_win);
    mov_res_win = SeeResiduals(mov_res_win, mc_win.^2);
    mov_res_win = SeeResiduals(mov_res_win, mc_win(:,1) .* mc_win(:,end));

    % 5. Bleaching + blood vessel
    mov_res_win = SeeResiduals(mov_res_win, bkg_win, 1);

    % 6. dF/F and slow-trend subtraction
    ref_im_sm   = imgaussfilt(VoltResult.ref_im, volt_smooth_sigma);
    dFF_win     = -imgaussfilt(mov_res_win, volt_smooth_sigma) ./ ref_im_sm;
    dFF_win_sub = movmean(dFF_win, 20, 3);

    % 7. Subtract pre-trigger baseline per pixel
    baseline     = mean(dFF_win(:,:, 1:N_pre), 3);
    dFF_win      = dFF_win     - baseline;
    dFF_win_sub  = dFF_win_sub - mean(dFF_win_sub(:,:, 1:N_pre), 3);

    %-- Accumulate into STA
    STA_volt     = STA_volt     + dFF_win;
    STA_volt_sub = STA_volt_sub + dFF_win_sub;
    nContrib     = nContrib     + 1;

    fprintf('  Event %d / %d done.\n', ei, nEvents);
end

%-- Normalise STA
STA_volt     = STA_volt     ./ reshape(nContrib, 1, 1, []);
STA_volt_sub = STA_volt_sub ./ reshape(nContrib, 1, 1, []);

%% ===================== WARP STA INTO GLU SPACE ==========================
% Transform one test frame to get output dimensions
test_warp = transformCamera_O2B(Device_Data, tformReg, STA_volt(:,:,1), GluResult.AvgGluImg);
[H_out, W_out] = size(test_warp);

%% ===================== FIND GLU FRAME NEAREST EACH VOLT FRAME ===========
% t_glu_show : time axis of mov_glu_dFF_masked
% For each tau offset, find the closest Glu frame to display alongside

t_trig_s     = glu_event_times_s(valid);   % trigger times in seconds
t_win_s_base = tau_axis / Fs_volt;         % relative time in seconds

%% ===================== VIDEO WRITER =====================================
v = VideoWriter(output_STA_filename, 'MPEG-4');
v.FrameRate = frame_rate;
v.Quality   = quality;
open(v);
fig = figure('Color','k','Units','normalized','Position',[0 0 1 0.45]);

%% ===================== RENDER LOOP ======================================
fprintf('Rendering STA movie...\n');

for fi = 1:nWin

    %-- Warp STA volt frames into Glu space
    sta_frame     = STA_volt(:,:,fi);
    sta_frame_sub = STA_volt_sub(:,:,fi);

    alignedSTA     = transformCamera_O2B(Device_Data, tformReg, sta_frame,     GluResult.AvgGluImg);
    alignedSTA_sub = transformCamera_O2B(Device_Data, tformReg, sta_frame_sub, GluResult.AvgGluImg);

    volt_rgb     = grs2rgb(alignedSTA,     cmap_volt, volt_dFF_range(1), volt_dFF_range(2));
    volt_sub_rgb = grs2rgb(alignedSTA_sub, cmap_volt, volt_dFF_range(1), volt_dFF_range(2));

    volt_rgb     = volt_rgb     .* DendriteMask;
    volt_sub_rgb = volt_sub_rgb .* DendriteMask;

    %-- Get corresponding Glu frame (use mean trigger time as reference)
    t_abs_s   = mean(glu_event_times_s(valid)) + t_win_s_base(fi);
    gi_frac   = interp1(t_glu_show, 1:numel(t_glu_show), t_abs_s, 'linear', 'extrap');
    gi_frac   = min(max(gi_frac, 1), size(mov_glu_dFF_masked, 3));
    gi_lo     = max(floor(gi_frac), 1);
    gi_hi     = min(ceil(gi_frac),  size(mov_glu_dFF_masked, 3));
    alpha_g   = gi_frac - gi_lo;
    glu_frame = (1-alpha_g) * mov_glu_dFF_masked(:,:,gi_lo) + alpha_g * mov_glu_dFF_masked(:,:,gi_hi);
    glu_rgb   = grs2rgb(glu_frame, cmap_glu, glu_dFF_range(1), glu_dFF_range(2));

    %-- Overlay
    overlay = glu_rgb + volt_sub_rgb;
    overlay = min(overlay, 1);

    %-- Render
    clf(fig);
    set(fig, 'Color','k');
    tl = tiledlayout(fig, 1, 3, 'TileSpacing','compact', 'Padding','compact');

    ax1 = nexttile(tl, 1);
    imshow(glu_rgb, 'Parent', ax1);
    title(ax1, sprintf('Glutamate  (%+.0f ms)', t_axis_ms(fi)), 'Color','w', 'FontSize',11);

    ax2 = nexttile(tl, 2);
    imshow(volt_rgb, 'Parent', ax2);
    title(ax2, sprintf('Volt STA  n=%d  (%+.0f ms)', nEvents, t_axis_ms(fi)), 'Color','w', 'FontSize',11);

    ax3 = nexttile(tl, 3);
    imshow(overlay, 'Parent', ax3);
    title(ax3, 'Overlay (G=Glu, R=Volt STA)', 'Color','w', 'FontSize',11);

    drawnow;
    writeVideo(v, getframe(fig));

    if mod(fi, 50) == 0
        fprintf('  Frame %d / %d  (%.1f%%)\n', fi, nWin, 100*fi/nWin);
    end
end

close(v);
close(fig);
fprintf('Done. STA movie saved to:\n  %s\n', output_STA_filename);
