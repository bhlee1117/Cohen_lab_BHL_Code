%% #########################################################################
%  GLUTAMATE-TRIGGERED AVERAGE VOLTAGE  —  ANALYSIS PIPELINE
%
%  Flow (each %% block is a runnable section):
%    0. Paths & user settings
%    1. Load data & camera geometry
%    2. Shared preprocessing (masks, footprints, spikes, traces, distances)
%    3. Trace analysis   — calculations
%    4. Trace analysis   — figures
%    5. Movie analysis   — calculations
%    6. Movie analysis   — figures & movie export
%    7. Interactive explorer
% #########################################################################

%% ===== 0. PATHS & USER SETTINGS =========================================
clear; clc;
cd '/Volumes/BHL18TB_D2/Arranged_Data/Prism_OptopatchResult';
[~, ~, raw] = xlsread( ...
    '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement.xlsx', ...
    'Sheet1', 'C5:AA31');
load('/Volumes/BHL18TB_D2/20260203_SD_V2+iGluSNFR4/25X_transformationMatrix.mat');  % -> tformReg
fpath         = raw(:,1)';
StructureData = raw(:,6);

% --- Session / acquisition ---
f            = 3;        % session index into fpath{}
framesPerSeg = 15000;    % frames per mc0x.bin file (15 s x 1000 Hz)
Fs_volt      = 1000;     % voltage camera frame rate (Hz)
pixelsize    = 0.4680;   % um per pixel

% --- Movie loading limits ---
V2moviemaxTime  = 15000;
GlumoviemaxTime = 5000;

% --- Spatial smoothing (Gaussian sigma in pixels) ---
volt_smooth_sigma = 1.0;
glu_smooth_sigma  = 1.5;
strImg_bin_volt   = [];  % optional structure mask (H_v x W_v, 0-1); [] to skip

% --- Analysis parameters ---
VoltLag    = 1000;          % voltage STA half-window for trace analysis (ms)
caxis      = [-0.005 0.01]; % color limits for voltage STA panels
GluPeakWin = 2;             % peri-spike half-window (glu frames) for sub-frame peak spline

%% ===== 1. LOAD DATA & CAMERA GEOMETRY ===================================
load(fullfile(fpath{f}, 'output_data.mat'));            % -> Device_Data
GluResult  = importdata(fullfile(fpath{f},"Glu_Result.mat"));
VoltResult = importdata(fullfile(fpath{f},"Volt_Result.mat"));

% ---- Correct time axes: t_ax is stored one frame longer than the traces ----
% Trim each time axis to the number of samples in its signal trace so that
% t_ax and the traces are indexed consistently everywhere downstream.
GluResult.t_ax  = GluResult.t_ax(1:size(GluResult.dFF_glu, 2));
VoltResult.t_ax = VoltResult.t_ax(1:size(VoltResult.normTraces, 2));

% Voltage camera geometry
nCol1 = double(Device_Data{3}.ROI(2));
nRow1 = double(Device_Data{3}.ROI(4));
exposuretime1 = 1 / Fs_volt;                   % s per frame (= 0.001 s)

% Glutamate camera geometry
nCol2 = double(Device_Data{4}.ROI(2));
nRow2 = double(Device_Data{4}.ROI(4));
exposuretime2 = Device_Data{4}.exposuretime;   % s per frame (slower camera)

%% ===== 2. SHARED PREPROCESSING ==========================================
% ---- Dendrite tracing -> masks (for geodesic distances & movie overlay) ----
StructureStack = read_tiff(StructureData{f});
sz_struct = size(StructureStack(:,:,1));
[H_out, W_out] = size(GluResult.AvgGluImg);

swcfiles = dir(fullfile(fpath{f}, 'Tracing*.swc'));
if isempty(swcfiles)
    swcfiles = dir(fullfile(fileparts(StructureData{f}), 'Tracing*.swc'));
end
dendritePoints = [];
for si = 1:length(swcfiles)
    swcname = fullfile(swcfiles(si).folder, swcfiles(si).name);
    tree = load_tree(swcname);
    dendritePoints{si} = [tree.X, tree.Y];
end
DendriteMask        = point2img(cell2mat(dendritePoints'),  3.5, sz_struct);
P_str = cell2mat(dendritePoints(1));          % [x y]
DendriteMask_N1_str = point2img(P_str , 3,  sz_struct);
DendriteMask_N1     = imwarp(DendriteMask_N1_str, VoltResult.tform_Str2Volt, ...
                             'OutputView', imref2d(size(VoltResult.ref_im)));

Rvolt = refFromROI(size(VoltResult.ref_im),   double(Device_Data{1,3}.ROI));
Rglu  = refFromROI(size(GluResult.AvgGluImg), double(Device_Data{1,4}.ROI));
[xv , yv ] = transformPointsForward(VoltResult.tform_Str2Volt, P_str(:,1), P_str(:,2));
[xvw, yvw] = intrinsicToWorld(Rvolt, xv, yv);
[xgw, ygw] = transformPointsInverse(tformReg, xvw, yvw);   % Volt world -> Glu world
[xg , yg ] = worldToIntrinsic(Rglu, xgw, ygw);
DendriteMask_N1_glu = point2img([xg yg], 3, size(GluResult.AvgGluImg));

% ---- Voltage reference image & footprints in glutamate-camera space ----
VoltRefIm     = transformCamera_O2B(Device_Data, tformReg, VoltResult.ref_im, GluResult.AvgGluImg);
VoltRefIm_rgb = grs2rgb(VoltRefIm, gray(256), 20, 200);
AlignedVolt_ftprnt = transformCamera_O2B(Device_Data, tformReg, VoltResult.ftprnt, GluResult.AvgGluImg);
nROI = size(VoltResult.ftprnt, 3);

% ---- Glutamate spikes & positions ----
GluSpike = find_spike_bh(zscore(GluResult.dFF_glu,0,2), 3, 1);
SumTglu  = sum(GluSpike, 1);
SumGlu   = sum(GluSpike, 2);
[~, synsort] = sort(SumGlu, 'ascend');
%show_traces_spikes(GluResult.dFF_glu(synsort,:),GluSpike(synsort,:),SumTglu)
Glu_coord = get_coord(GluResult.S_glu);

% ---- Normalized voltage traces (per-ROI dF/F) ----
F0_volt        = tovec(VoltResult.ref_im)' * tovec(VoltResult.ftprnt);
VoltTrace_Norm = VoltResult.normTraces ./ F0_volt';

% ---- Glutamate frame timing ----
GluTimeGap = min(diff(GluResult.t_ax));
GluLag     = ceil(VoltLag / (GluTimeGap*1000));

% ---- Dendrite mask for geodesic distances ----
% Keep a full-mask copy (SkelDend) for the per-synapse geodesic distances that
% are computed in section 4 (once Syn2show is chosen); then remove somatic ROIs
% from DendriteMask_N1 for all downstream masking/overlays.

if isfile(fullfile(fpath{f},'DendriteMaskN1.mat'))
    load(fullfile(fpath{f},'DendriteMaskN1.mat'))
SkelDend = DendriteMask_N1_glu;
else
[~, rmvROI] = get_ROI(DendriteMask_N1_glu);
DendriteMask_N1_glu= DendriteMask_N1_glu.* (max(rmvROI,[],3)==0);
SkelDend = DendriteMask_N1_glu;
DendriteMask_N1= transformCamera(Device_Data, tformReg, VoltResult.ref_im, DendriteMask_N1_glu);
save(fullfile(fpath{f},'DendriteMaskN1.mat'),'DendriteMask_N1','DendriteMask_N1_glu','-v7.3');
end

%% ===== 2.5. Movie PREPROCESSING ==========================================
% ---- (a) Glutamate movie frame timing (rebuild cam2 trigger in seconds) ----
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

t_glu = GluResult.t_ax(1:size(GluResult.dFF_glu,2));   % glu frame times (s)

if cam2_vsyn(end) < GlumoviemaxTime
    GluemovTimesegments = [cam2_vsyn(mod488)+2, cam2_vsyn(end)];
else
    GluemovTimesegments = [(cam2_vsyn(mod488)+2) : GlumoviemaxTime : cam2_vsyn(end)];
    GluemovTimesegments(end+1) = cam2_vsyn(end);
end

% ---- (b) Load & process glutamate movie -> per-ROI dF/F image (dFF_roi) ----
n_kernel = 1.5;   % median filter kernel size (pixels)
dil_rad  = 2;     % dilation radius (pixels)

fprintf('Loading glutamate movie segments...\n');
mov_glu = readBinMov_BHL_multiple(fpath{f}, 4, 1:length(t_glu), 5000, 'mc2');
szGlu   = size(mov_glu);

mov_glu2 = SeeResiduals(mov_glu , GluResult.mc);
mov_glu2 = SeeResiduals(mov_glu2, GluResult.mc.^2);
mov_glu2 = SeeResiduals(mov_glu2, GluResult.mc(:,1).*GluResult.mc(:,end));

[V, D, u_sub] = get_eigvector(tovec(imgaussfilt(mov_glu2,glu_smooth_sigma))', 40);
u_sub = reshape(u_sub, szGlu(1), szGlu(2), []);
figure(2); clf; imshow2_patch(u_sub); drawnow;
n = input("PCs to regress out\n");
mov_glu2 = SeeResiduals(mov_glu2, V(:,n));

N_roi = size(GluResult.S_glu, 3);

% Normalize by mean of raw mov_glu within each ROI
F0_roi = tovec(mean(mov_glu,3))' * tovec(GluResult.S_glu) ./ squeeze(sum(GluResult.S_glu,[1 2]))';
roi_mask     = max(GluResult.S_glu>0, [], 3);
roi_mask_dil = imdilate(roi_mask, strel('disk', dil_rad));
mov_roi      = maskBinary(mov_glu2, roi_mask_dil==0, NaN);
mov_roi      = imgaussfiltnan(mov_roi, n_kernel);
F0_roi_img   = zeros(size(mov_glu2(:,:,1)), 'double');
for i = 1:N_roi
    roi_mask = GluResult.S_glu(:,:,i) > 0;
    roi_mask = imdilate(roi_mask, strel('disk', dil_rad));
    F0_roi_img_tmp = maskBinary(zeros(size(mov_glu2(:,:,1)),'double'), roi_mask, F0_roi(i));
    F0_roi_img     = max(cat(3, F0_roi_img, F0_roi_img_tmp), [], 3);
end
dFF_roi = (mov_roi) ./ F0_roi_img .* max(GluResult.S_glu>0, [], 3);
clear mov_roi mov_glu u_sub

%% ===== 3. TRACE ANALYSIS — CALCULATIONS =================================
% ---- (a) Peri-spike matrices for EVERY glutamate ROI ----
% Computed for all synapses and stored in cell arrays; the synapse to display
% (Syn2show) is selected later, in section 4. Trigger times are refined to
% sub-frame precision (see gluPeaks2VoltFrames at the end of this file): the
% glutamate signal is sampled slower than voltage, so instead of snapping each
% spike to the nearest voltage frame we spline-interpolate the peri-spike
% glutamate signal onto the voltage time axis and trigger on the interpolated
% peak, so the voltage average is better aligned.
%
% NB: gspMat_voltage_all stores a +/-VoltLag window per event for every synapse,
%     which can be memory-heavy; reduce VoltLag or the stored window if needed.
nGluSyn            = size(GluSpike,1);
gspMat_all         = cell(nGluSyn,1);   % glu-spike raster STA  [nGluROI x nEv x (2*GluLag+1)]
gspMat_voltage_all = cell(nGluSyn,1);   % voltage STA per event [nROI    x nEv x (2*VoltLag+1)]
for s = 1:nGluSyn
    sf = find(GluSpike(s,:) > 0);       % glutamate spike frames for this synapse
    gsptime_volt = gluPeaks2VoltFrames(sf, GluResult.dFF_glu(s,:), ...
                                       GluResult.t_ax, VoltResult.t_ax, GluPeakWin);
    [~, gmv, validglu] = get_STA(VoltTrace_Norm, gsptime_volt, VoltLag, VoltLag);
    validglu = ismember(gsptime_volt, validglu);
    gspike_tmp=find(GluSpike(s,:)>0);
    [~, gm]            = get_STA(GluSpike, gspike_tmp(validglu), GluLag, GluLag);
    gspMat_all{s}         = gm;
    gspMat_voltage_all{s} = gmv;
end
tax_glu_sub = [-GluLag:GluLag]*GluTimeGap;

% ---- (b) Spatiotemporal clustering & subthreshold vs cluster size ----
tau = 1; dist_threshold = 30; %um
GluMatched2Neuron = matchSglu2SWC(dendritePoints, Glu_coord, 'MaxDist', 5);
Syn_N1        = GluMatched2Neuron.neuronID==1;
GluSpike_sub  = GluSpike(Syn_N1,:);
Glu_coord_sub = Glu_coord(Syn_N1,:);
nSyn          = sum(Syn_N1);

% Pairwise geodesic distances between N1 synapses (precomputed & loaded).
load('/Volumes/BHL18TB_D2/20260203_SD_V2+iGluSNFR4/170544BHLm302_V_Glu_25X_3min_N2/DistanceMatrix_GluSignal_N1.mat')  % -> SynDistance_sub
% SynDistance_sub=[];
% for i=1:nSyn
%     for j=i:nSyn
%         SynDistance_sub(i,j)=geodesic_distance(SkelDend,Glu_coord_sub(i,:),Glu_coord_sub(j,:));
%     end
%     if mod(i,20)==0, fprintf('  %d / %d (%.1f%%)\n', i, nSyn, 100*i/nSyn); end
% end
% SynDistance_sub=max(cat(3,SynDistance_sub,SynDistance_sub'),[],3);
% save('/Volumes/BHL18TB_D2/20260203_SD_V2+iGluSNFR4/170544BHLm302_V_Glu_25X_3min_N2/DistanceMatrix_GluSignal_N1.mat','SynDistance_sub','-v7.3')  % -> SynDistance_sub

GluEvents_cluster = find_spatiotemporal_events(GluSpike_sub, SynDistance_sub, Glu_coord_sub, ...
                                               'MaxDist', 10/pixelsize, 'Tolerance', 1);
EventCentroid          = vertcat(GluEvents_cluster.centroid);
CentroidFtprnt_distMat = distance_mat(EventCentroid, get_coord(AlignedVolt_ftprnt));
[Cent2Ftprnt_dist, minarg] = min(CentroidFtprnt_distMat, [], 2);
validEvent         = Cent2Ftprnt_dist*pixelsize < dist_threshold;
GluEvents_cluster2 = GluEvents_cluster(find(validEvent));
minarg             = minarg(validEvent);
GluFtprnt_distMat  = distance_mat(Glu_coord(Syn_N1,:), get_coord(AlignedVolt_ftprnt));
[Glu2Ftprnt_dist, minarg2] = min(GluFtprnt_distMat, [], 2);

[eventinVolt] = match_nearest(GluResult.t_ax(vertcat(GluEvents_cluster(find(validEvent)).time)), VoltResult.t_ax);
[~, eventMat_volt, eventmatchedVoltageTime] = get_STA(VoltTrace_Norm, eventinVolt, VoltLag, VoltLag);
eventMat_volt      = eventMat_volt - median(eventMat_volt, 3);
isVoltValid        = ismember(eventinVolt, eventmatchedVoltageTime);
GluEvents_voltValid = GluEvents_cluster2(isVoltValid);

EventclusterSize = vertcat(GluEvents_voltValid.size);
minarg           = minarg(isVoltValid);
EventVoltage     = [];
for ev = 1:size(eventMat_volt,2)
    EventVoltage(ev,:) = eventMat_volt(minarg(ev),ev,:);
end
EventVoltageClustered = [];
for clsz = [min(EventclusterSize):max(EventclusterSize)]
    EventVoltageClustered{clsz} = mean(EventVoltage(EventclusterSize==clsz, VoltLag+[-30:30]), 2, 'omitnan');
end

% Depolarization variation across synapses at isolated events
roi_event_idx  = cell(nSyn, 1);    % event indices for each ROI
roi_event_size = cell(nSyn, 1);    % event sizes   for each ROI
for n = 1:nSyn
    roi_event_idx{n}  = find(arrayfun(@(e) ismember(n, e.rois), GluEvents_voltValid));
    roi_event_size{n} = [GluEvents_voltValid(roi_event_idx{n}).size];
end
clustersize2show   = [1];
median_clustersize = cellfun(@(x) prctile(x,90), roi_event_size);
idx2show           = cellfun(@(x,y) y(ismember(x,clustersize2show)), roi_event_size, roi_event_idx, 'UniformOutput', false);
IsolatedSubthreshold = [];
for n = 1:nSyn
    subthvec = mean(mean(eventMat_volt(minarg2(n), idx2show{n}, VoltLag+[-30:30]), [3]), 2);
    IsolatedSubthreshold = [IsolatedSubthreshold; [n median_clustersize(n) subthvec]];
end

%% ===== 3.5. INTERACTIVE GLU-TRIGGERED EXPLORER ============================
% Click a glutamate synapse to explore its glutamate-triggered averages and
% peri-spike matrices. See interactive_glutamateTA_viewer.m for details.
interactive_glutamateTA_viewer(GluResult, VoltResult, VoltTrace_Norm, ...
    GluSpike, AlignedVolt_ftprnt, VoltRefIm, DendriteMask_N1, Device_Data, tformReg, ...
    'PixelSize', pixelsize, 'VoltLag', 100, 'GluWin_ms', 200, 'InitROI', 7);  % initial ROI only

%% ===== 3.5b. NOMINATED ROIs — EVENT TABLE + GLU-TRIGGERED AVERAGE VOLTAGE =====
% For a set of nominated (good) glutamate ROIs, this block:
%   (1) builds a per-event table (glu ROI, position, glu time, spline-interpolated
%       glu time, nearest voltage time, # coactive glutamate in neighbours, and
%       the voltage subthreshold amplitude), and
%   (2) extracts the peri-event voltage movie for each event, concatenates and
%       saves it to .bin files (same pattern as Analysis_20260518_SS_CS_...),
%       and computes the glu-triggered average voltage (trace kymo + movie).
% NB: this is a self-contained analysis; it does not feed section 5 (which is a
%     separate single-synapse Syn2look movie).

% ===== USER: nominate good glutamate ROIs =====
NominatedROIs = [7 8 10 64 99 266 109 259 183 241 167];            % <-- edit: indices into GluResult.S_glu / GluSpike

% ---- Parameters ----
SplineUp    = 20;               % spline upsampling for the sub-frame glu peak time
coactRadius = 50;               % um, neighbourhood radius for the coactive-glutamate count
coactTol    = 1;                % glu-frame tolerance for coactivity (+/- frames)
subthWin_ms = 20;               % +/- window (ms) to average the voltage subthreshold at the event
movHalf     = 100;             % +/- window (ms) for the peri-event voltage movie / STA
tau_axis    = -movHalf : movHalf;
nWin        = numel(tau_axis);
centerIdx   = movHalf + 1;      % index of t = 0 in the window
nGluFrame   = size(GluSpike, 2);
nVoltFrame  = numel(VoltResult.t_ax);

% ---- Nearest voltage ROI for every glutamate ROI (euclidean, glu space) ----
[~, nearestVoltROI_all] = min(distance_mat(Glu_coord, get_coord(AlignedVolt_ftprnt)), [], 2);

% ---- Geodesic glu-glu distances from nominated ROIs to all ROIs (coactivity) ----
% Same convention as SynDistance_sub (geodesic on SkelDend with glu-space coords).
GluDistFromNom = NaN(numel(NominatedROIs), size(GluSpike,1));
for a = 1:numel(NominatedROIs)
    fprintf('Calculating geodesic distance of ROI# %2.0f \n',NominatedROIs(a));
    for b = find(Syn_N1)'
        GluDistFromNom(a,b) = geodesic_distance(SkelDend, Glu_coord(NominatedROIs(a),:), Glu_coord(b,:));        
    end
end

% ---- Build the per-event list ----
evGluROI=[]; evX=[]; evY=[]; evGluTime=[]; evInterpTime=[]; evVoltTime=[]; 
evVoltFrame=[]; evNcoact=[]; evGluAmp=[]; evGluFrame=[];
for a = 1:numel(NominatedROIs)
    s     = NominatedROIs(a);
    sf    = find(GluSpike(s,:) > 0);                         % glutamate spike frames
    neigh = setdiff(find(GluDistFromNom(a,:)*pixelsize < coactRadius), s);   % neighbours within radius
    for t = sf
        % sub-frame interpolated glutamate peak time (spline on a fine grid)
        w  = max(1,t-GluPeakWin) : min(nGluFrame, t+GluPeakWin);
        tl = GluResult.t_ax(w);  yl = GluResult.dFF_glu(s,w);
        if numel(tl) >= 3
            tq = linspace(tl(1), tl(end), (numel(tl)-1)*SplineUp + 1);
            yq = interp1(tl, yl, tq, 'spline');
            [~, im] = max(yq);  t_interp = tq(im);
        else
            t_interp = GluResult.t_ax(t);
        end
        vfr  = match_nearest(t_interp, VoltResult.t_ax);     % nearest voltage frame to interp. peak
        twin = max(1,t-coactTol) : min(nGluFrame, t+coactTol);
        nco  = sum(any(GluSpike(neigh, twin) > 0, 2));       % # coactive neighbour glutamate

        evGluROI(end+1,1)     = s;
        evX(end+1,1)          = Glu_coord(s,1);
        evY(end+1,1)          = Glu_coord(s,2);
        evGluTime(end+1,1)    = GluResult.t_ax(t);
        evGluFrame(end+1,1) = t;
        evInterpTime(end+1,1) = t_interp;
        evVoltFrame(end+1,1)  = vfr;
        evVoltTime(end+1,1)   = VoltResult.t_ax(vfr);
        evNcoact(end+1,1)     = nco;
        evGluAmp(end+1,1) = GluResult.dFF_glu(s,t);
    end
end

% ---- Keep events whose full movie window is inside the recording ----
validEv = (evVoltFrame + tau_axis(1) >= 1) & (evVoltFrame + tau_axis(end) <= nVoltFrame);
evGluROI=evGluROI(validEv); evX=evX(validEv); evY=evY(validEv);
evGluTime=evGluTime(validEv); evInterpTime=evInterpTime(validEv);
evVoltTime=evVoltTime(validEv); evVoltFrame=evVoltFrame(validEv); evNcoact=evNcoact(validEv);
nEv = numel(evVoltFrame);
fprintf('3.5b: %d nominated events across %d ROI(s).\n', nEv, numel(NominatedROIs));

% ---- Glu-triggered average voltage (trace) + per-event subthreshold amplitude ----
[~, MatV] = get_STA(VoltTrace_Norm, evVoltFrame', movHalf, movHalf);   % nROI x nEv x nWin
MatV = MatV - median(MatV, 3);                                          % baseline per ROI/event
MatV = movmean(MatV,20,3);
GluShift=NaN(length(evGluROI),1);
for n=NominatedROIs
    evn =find(evGluROI==n);
    nv = nearestVoltROI_all(n);
    GTA =squeeze(mean(MatV(nv, evn , :),2,'omitnan'));
    [~, Gsh]=max(GTA(movHalf+[-round(GluTimeGap*1000):round(GluTimeGap*1000)]));
    GluShift(evn)= -round(GluTimeGap*1000)+Gsh;
end

evSubthr = zeros(nEv,1);
for e = 1:nEv
    nv = nearestVoltROI_all(evGluROI(e));                               % nearest voltage ROI
    evSubthr(e) = mean(MatV(nv, e, movHalf + GluShift(e)+1 + (-subthWin_ms:subthWin_ms)), 3);
end
GluTrigVoltSTA = squeeze(mean(MatV, 2, 'omitnan'));                     % nROI x nWin (trace kymo)

% ---- Assemble & save the event table ----
EventTable = table(evGluROI, evX, evY, evGluTime, evInterpTime, evVoltTime, evNcoact, evSubthr, evGluAmp, evVoltFrame, GluShift, evGluFrame, ...
    'VariableNames', {'GluROI','GluX','GluY','GluTime_s','InterpGluTime_s', ...
                      'NearestVoltTime_s','nCoactive','VoltSubthr','GluAmp','VoltFrame','SubFrameGluShift','GluFrame'});
disp(EventTable);
save(fullfile(fpath{f}, 'GluEventTable_nominated.mat'), ...
     'EventTable', 'GluTrigVoltSTA', 'tau_axis', 'NominatedROIs', '-v7.3');

%% ===== 3.5b (movie). PERI-EVENT VOLTAGE MOVIE — CONCATENATE, SAVE, AVERAGE =====
% ---- Voltage movie regressors (bleaching profile + Fstd image) ----
fprintf('Preparing voltage movie regressors...\n');
volt_smooth_sigma = 1;
t2bleachfit = 15001:25001;
mov_ref_seg = readBinMov_BHL_multiple(fpath{f}, 3, t2bleachfit, framesPerSeg, 'mc');
meanF_ref   = squeeze(mean(mov_ref_seg, [1 2]));
bleach_ref  = expfitDM_2(t2bleachfit', meanF_ref, (1:nVoltFrame)', 10000);
bleach_ref  = bleach_ref / mean(bleach_ref);                           % unit-mean bleach profile
% mc_win = VoltResult.mc(t2bleachfit, :);
% mov_ref_seg = SeeResiduals(mov_ref_seg, bleach_ref(t2bleachfit));
% mov_ref_seg = SeeResiduals(mov_ref_seg, mc_win);
% mov_ref_seg = SeeResiduals(mov_ref_seg, mc_win.^2);
% mov_ref_seg = SeeResiduals(mov_ref_seg, mc_win(:,1) .* mc_win(:,end));
% mov_ref_seg = SeeResiduals(mov_ref_seg, VoltResult.bvTrace(:, t2bleachfit), 1);
% Fstd_img = get_F0img_PCA(imgaussfilt(mov_ref_seg(:,:,1:5:10000), volt_smooth_sigma));
% clear mov_ref_seg
sz_volt = size(VoltResult.ref_im);
nPix    = sz_volt(1) * sz_volt(2);

nSeg           = ceil(nVoltFrame / framesPerSeg);
movfile_prefix = fullfile(fpath{f}, 'GluTrigger_VoltMovie_');
size_limit     = 2.5e9;                              % bytes; flush to a new .bin beyond this
MovieInfo = struct('voltFrame',{{}}, 'gluROI',{{}}, 'StackedMovieN',0);
Mov_PeakTA = [];  c = 1;  MovieInfo.voltFrame{c} = [];  MovieInfo.gluROI{c} = [];
GluTrigVoltMovie_sum = zeros(nPix, nWin);            % running sum of raw-dF windows

for j = 1:nSeg
    core_s = (j-1)*framesPerSeg + 1;                 % this chunk's frames (no overlap)
    core_e = min(j*framesPerSeg, nVoltFrame);
    evIdx  = find(evVoltFrame >= core_s & evVoltFrame <= core_e);   % events centered in this chunk
    if isempty(evIdx)
        fprintf('  chunk %d/%d: no events, skipping.\n', j, nSeg);
        continue;
    end
    pad_s    = max(1, core_s - movHalf);             % pad so peri-event windows near edges are covered
    pad_e    = min(nVoltFrame, core_e + movHalf);
    padRange = pad_s:pad_e;
    fprintf('  chunk %d/%d: %d events, loading frames %d-%d ...\n', j, nSeg, numel(evIdx), pad_s, pad_e);

    % 1) load the chunk (padded)
    mov_seg = readBinMov_BHL_multiple(fpath{f}, 3, padRange, framesPerSeg, 'mc');

    % 2) regress the whole chunk (same operations as the per-event version)
    mov_seg = mov_seg - mean(mov_seg, 3);
    bkg     = [bleach_ref(padRange);];
    mc_seg  = VoltResult.mc(padRange, :);
    mov_seg = SeeResiduals(mov_seg, mc_seg);
    mov_seg = SeeResiduals(mov_seg, mc_seg.^2);
    mov_seg = SeeResiduals(mov_seg, mc_seg(:,1) .* mc_seg(:,end));
    mov_seg = SeeResiduals(mov_seg, bkg, 1);
    mov_seg = tovec(mov_seg);                         % pixels x T(chunk)

    % 3) take the frames we need (peri-glutamate windows) in chunk-local indices
    VoltFrame2add_shifted=EventTable.VoltFrame(evIdx) + EventTable.SubFrameGluShift(evIdx);
    localIdx =  VoltFrame2add_shifted - pad_s + 1;
    [~, AddMov, kept] = get_STA(mov_seg, localIdx', movHalf, movHalf);   % pixels x nEv x nWin
    addedInd=ismember(localIdx,kept);
    if numel(kept) ~= numel(localIdx)
        warning('chunk %d: %d / %d windows dropped at edges', j, numel(localIdx)-numel(kept), numel(localIdx));
    end
    clear mov_seg

    %GluTrigVoltMovie_sum = GluTrigVoltMovie_sum + squeeze(sum(AddMov, 2));   % accumulate for the average
    Mov_PeakTA = cat(3, Mov_PeakTA, permute(AddMov, [1 3 2]));               % pixels x time x events
    MovieInfo.voltFrame{c} = [MovieInfo.voltFrame{c} VoltFrame2add_shifted(addedInd)'];
    MovieInfo.gluROI{c}    = [MovieInfo.gluROI{c}    EventTable.GluROI(evIdx(addedInd))'];

    % 4) flush to a .bin file once the stack gets big
    Movinfo = whos('Mov_PeakTA');
    if Movinfo.bytes > size_limit
        MovtoWrite = vm(double(Mov_PeakTA) + 10000);
        MovtoWrite.transpose.savebin([movfile_prefix num2str(c) '.bin']);
        MovieInfo.StackedMovieN = MovieInfo.StackedMovieN + size(Mov_PeakTA, 3);
        c = c + 1;  MovieInfo.voltFrame{c} = [];  MovieInfo.gluROI{c} = [];
        clear Mov_PeakTA MovtoWrite;  Mov_PeakTA = [];
        fprintf('  saved movie chunk file %d\n', c-1);
    end
end
% flush the remaining stacked events
if ~isempty(Mov_PeakTA)
    MovtoWrite = vm(double(Mov_PeakTA) + 10000);
    MovtoWrite.transpose.savebin([movfile_prefix num2str(c) '.bin']);
    MovieInfo.StackedMovieN = MovieInfo.StackedMovieN + size(Mov_PeakTA, 3);
    clear Mov_PeakTA MovtoWrite;  Mov_PeakTA = [];
    fprintf('  saved movie chunk file %d\n', c);
end
MovieInfo.voltFrame(cellfun(@isempty, MovieInfo.voltFrame)) = [];
MovieInfo.gluROI(cellfun(@isempty, MovieInfo.gluROI))       = [];
save(fullfile(fpath{f}, 'GluTrigVoltMovieInfo.mat'), 'MovieInfo', 'tau_axis', 'movfile_prefix', 'sz_volt', '-v7.3');
%% ===== 3.5c (movie). Load glu-triggered movie and average
load(fullfile(fpath{f}, 'GluEventTable_nominated.mat'))
load(fullfile(fpath{f}, 'GluTrigVoltMovieInfo.mat'))
g=1; GTAmovie=[]; GTAGluTrace=[]; GTAvoltTrace=[]; GTAGluSpike=[]; GTAgmovie=[]; distSyn2Volt=[];
for n=unique(EventTable.GluROI)'
    for nv=1:size(AlignedVolt_ftprnt,3)
    [distSyn2Volt(nv,g), p]=geodesic_distance(SkelDend,Glu_coord(n,:),get_coord(AlignedVolt_ftprnt(:,:,nv)));
    end
    nEV=find(EventTable.GluROI==n);
    framelist=EventTable.VoltFrame(nEV) + EventTable.SubFrameGluShift(nEV);
    Movmatch=cellfun(@(x) find(ismember(x,framelist)),MovieInfo.voltFrame,'UniformOutput',0);
    AlignMov=zeros(sz_volt(2)*sz_volt(1),length(tau_axis));
    for c=1:length(Movmatch)
        if ~isempty(Movmatch{c})
            movfile2read=[movfile_prefix num2str(c) '.bin'];
            try
            Movreadsub=(double(readBinMov_times(movfile2read,sz_volt(2)*sz_volt(1),length(tau_axis),Movmatch{c}))-10000);
            disp('readBinMov_times failed, reading one-by-one');
            catch
                Movreadsub=[];
                for cc=1:length(Movmatch{c})
            Movreadsub(:,:,cc)=(double(readBinMov_times(movfile2read,sz_volt(2)*sz_volt(1),length(tau_axis),Movmatch{c}(cc)))-10000);    
                end
            end
            AlignMov=AlignMov+sum(Movreadsub,3);
        end
    end
    GTAmovie(:,:,:,g)=reshape(AlignMov,sz_volt(1),sz_volt(2),[])./length(framelist);
[GTAvoltTrace(:,:,g)]=get_STA(VoltResult.normTraces,framelist',-tau_axis(1),tau_axis(end));
[GTAGluTrace(:,:,g)]=get_STA(GluResult.dFF_glu,EventTable.GluFrame(nEV)',10,10);
[GTAGluSpike(:,:,g)]=get_STA(GluSpike,EventTable.GluFrame(nEV)',2,2);
[GTAgmovie(:,:,:,g)]=toimg(get_STA(tovec(dFF_roi),EventTable.GluFrame(nEV)',15,15),szGlu(1),szGlu(2));
fprintf('GTA movie of ROI # %2.0f has loaded\n',n);
g=g+1;
end
GTAmovie=GTAmovie-mean(GTAmovie(:,:,[1:(-tau_axis(1))/2 end+tau_axis(1)/2:end],:),3);

%% Show Voltage + Glutamate
Show_GluROI=[7 10 21 36 64 74 79 91 302 631]; dffscale=0.7; cmapscatter=nebula(length(Show_GluROI));
g2show=3; VoltROI2see=12; interDD=[]; rmvROI=[23 24 25];
if ~isfield(VoltResult,'interDendDist')
    SkelDend_V=transformCamera(Device_Data, tformReg, VoltResult.ref_im, SkelDend);
    Vcoord=get_coord(VoltResult.ftprnt);
    for ii=VoltROI2see%1:size(VoltResult.ftprnt,3)
        for jj=1:size(VoltResult.ftprnt,3)
            interDD(ii,jj)=geodesic_distance(SkelDend_V,Vcoord(ii,:),Vcoord(jj,:));
        end
    end
end
if ~isfile(fullfile(fpath{f},'GTAvoltagefitResult.mat'))
    for gg=[1:11]
        fprintf('Loading nominated syn. # %2.0f \n',gg)
        [dax, dsort]=sort(distSyn2Volt(:,gg)*pixelsize);
        VoltAmp=mean(GTAvoltTrace(:,-tau_axis(1)+[-10:10],gg),2)./F0;
        VoltAmp2show{g}=[dax VoltAmp(dsort)];
        fitResult{g}=interactive_expfit(dax,VoltAmp(dsort));
        if fitResult{g}.aborted
            break;
        end
        g=g+1;
    end
    save(fullfile(fpath{f},'GTAvoltagefitResult.mat'),'fitResult','VoltAmp2show','-v7.3');
    disp('GTA voltage fit result has saved.')
else
    load(fullfile(fpath{f},'GTAvoltagefitResult.mat'));
end

noi=setdiff([1:size(VoltResult.ftprnt,3)],rmvROI);
F0=tovec(VoltResult.ref_im)'*tovec(VoltResult.ftprnt);
dFFGTA2F0=GTAvoltTrace(noi,:,g2show)./F0(noi)';
[dax, dsort]=sort(interDD(VoltROI2see,noi));

figure(20); clf;
tiledlayout(3,2,'TileSpacing','tight','Padding','compact');
nexttile(1,[1 1]);
imshow2(VoltRefIm,[20 140]);
xlim([250 950]); ylim([20 380]);
title('Voltron2-JF608')
nexttile(3,[1 1]);
imshow2(GluResult.ave_im2,[30 250]);
xlim([250 950]); ylim([20 380]);
title('iGluSnFR4f');
drawScaleBar(100/pixelsize,'horizontal','color','w','linewidth',3);
nexttile(5,[1 1]);
imshow2(VoltRefIm,[20 140]); hold all;
scatter(Glu_coord(:,1),Glu_coord(:,2),8,[0.8 0.6 0.6],'filled'); hold all;
scatter(Glu_coord(Show_GluROI,1),Glu_coord(Show_GluROI,2),30,cmapscatter,'filled'); 
xlim([250 950]); ylim([20 380]);
title('iGluSnFR4f ROIs')

nexttile(2,[3 1]);
l=plot(GluResult.t_ax,GluResult.dFF_glu(Show_GluROI,:)'-[1:length(Show_GluROI)]*dffscale);
arrayfun(@(l,c) set(l,'Color',c{:}),l,num2cell(cmapscatter,2));
drawScaleBar(60,'horizontal','color','k','linewidth',3,'position',[182 -7.5]);
drawScaleBar(1,'vertical','color','k','linewidth',3,'position',[182 -7.5]);
axis tight off;
set_font('Arial'); set_fontsize(14);
set_figsize(223,200);

figure(21); clf;
tiledlayout(2,2,'TileSpacing','tight','Padding','tight');
ax1=nexttile(1,[1 2]);
DendriteMask_V = imwarp(DendriteMask, VoltResult.tform_Str2Volt, 'OutputView', imref2d(size(VoltResult.ref_im)));
gluTAftprnt=movmean(-GTAmovie(:,:,:,g2show),10,3);
gluTAftprnt=max(imgaussfilt(gluTAftprnt(:,:,-tau_axis(1)+[-10:10]),1),[],3);
gluTAftprnt=grs2rgb(gluTAftprnt,turbo(256),0.3,3.5);
gluTAftprnt=gluTAftprnt.*mat2gray(VoltResult.ref_im)*1.5+mat2gray(VoltResult.ref_im)*0.5;
gluTAftprnt=transformCamera_O2B(Device_Data,tformReg,gluTAftprnt,GluResult.AvgGluImg);
imshow2(gluTAftprnt,[]); 
cb=colorbar; colormap(ax1,'turbo');
cb.Ticks=[0 1]; cb.TickLabels=[0.3 3.5]; cb.Label.String='-∆F';
xlim([250 950]); ylim([130 300]);
drawScaleBar(100/pixelsize,'horizontal','color','w');
title('Glutamate-triggered average voltage');

ax2=nexttile(3,[1 2]);
gluTAGftprnt=GTAgmovie(:,:,:,g2show);
gluTAGftprnt=max(gluTAGftprnt,[],3);
gluTAGftprnt_colored=grs2rgb(gluTAGftprnt,turbo(256),0.04,0.2);
gluTAGftprnt_colored=gluTAGftprnt_colored.*(gluTAGftprnt>0)+mat2gray(VoltRefIm)*0.5;
imshow2(gluTAGftprnt_colored,[]);
xlim([250 950]); ylim([130 300]);
cb=colorbar; colormap(ax2,'turbo');
cb.Ticks=[0 1]; cb.TickLabels=[4 20]; cb.Label.String='∆F/F (%)';
drawScaleBar(100/pixelsize,'horizontal','color','w');
title('Glutamate-triggered average glutamate');
set_font('Arial'); set_fontsize(18);
set_figsize(300,165);
%%
figure(22); clf;
tiledlayout(1,3,'TileSpacing','compact','Padding','tight');
ax3=nexttile(1,[1 1]);
imagesc([-100:100],[1:length(dsort)],dFFGTA2F0(dsort,:)*100,[-0.2 1.2]); 
set_kymoYtick(dax*pixelsize);
cb=colorbar; cb.Label.String='∆F/F (%)'; colormap(ax3,'turbo');
xlabel('Peri-glutamate time (ms)'); ylabel(sprintf('Distance from triggering \n iGluSnFR4f ROI (µm)')); xlabel('Time (ms)');
title('Glutamate-triggered average voltage')

nexttile(2,[1 1]);
peri_t_glu=[-10:10]*max(diff(GluResult.t_ax))*1000;
plot(peri_t_glu,GTAGluTrace(10,:,3)*100,'k','linewidth',1.5); box off;
xlim([-400 400]);
xlabel('Time (ms)'); ylabel('∆F/F (%)');
title('Glutamate-triggered average glutamate')

nexttile(3,[1 1]); hold all;
decayshowg=setdiff([1:11],[1 2 5 9 11]);
tconstGlu=cellfun(@(x) x.b,fitResult(decayshowg));
normVoltAmp2show=cellfun(@(x) [x(:,1) x(:,2)./max(x(:,2))],VoltAmp2show,'UniformOutput',0);
normVoltAmp2show2=[];
g=1;
for gg=decayshowg
    v2show=normVoltAmp2show{gg};
    v2show(v2show(:,2)<-0.2,:)=NaN;
    plot(v2show(:,1),v2show(:,2),'color',[0.6 0.6 0.6],'linewidth',1.5);
    normVoltAmp2show2{g}=[v2show];
    g=g+1;
end
binResult=binning_data(normVoltAmp2show2,[0:30:400]-15);
errorbar_shade(binResult.centers,binResult.mean,binResult.sem,[1 0 0]);
xlim([0 350]); ylim([-0.2 1.2]);
title(sprintf('L_c: %2.0f ± %2.0f µm', mean(tconstGlu),std(tconstGlu)))
xlabel(sprintf('Distance from triggering \n iGluSnFR4f ROI (µm)')); ylabel('Normalized voltage');

set_font('Arial'); set_fontsize(20);
set_figsize(400,115);

%% (Optional) SAVE as movie
g=3;
mov2write=-movmean(imgaussfilt(GTAmovie(:,:,:,g),1),10,3);
mov2write_rgb=[];
for i=1:size(mov2write,3)
mov2write_rgb(:,:,:,i)=grs2rgb(mov2write(:,:,i),turbo(256),0.3,3.5);
end
mov2write_rgb=mov2write_rgb.*mat2gray(VoltResult.ref_im)*1.5+mat2gray(VoltResult.ref_im)*0.5;
opt=[];
opt.timeVec=tau_axis;
opt.frameRate=15;
opt.colorbarLabel='-∆F';
opt.scaleBar=100/pixelsize;
opt.scaleBarText='100 µm';
opt.fontName='Arial';
opt.range=[0.3 3.5];
writeMov4d(fullfile(fpath{f},'GTA_movieV_Syn#10'),mov2write_rgb,opt)

mov2writeG=GTAgmovie(:,:,:,g);
mov2writeG_rgb=[];
for i=1:size(mov2writeG,3)
    glutmp=transformCamera(Device_Data,tformReg,VoltResult.ref_im,mov2writeG(:,:,i));
    glutmp(isnan(glutmp))=0;
mov2writeG_rgb(:,:,:,i)=grs2rgb(glutmp,turbo(256),0.01,0.1).*(glutmp>0);
end
mov2writeG_rgb=mov2writeG_rgb+mat2gray(VoltResult.ref_im)*0.5;
opt.timeVec=[-15:15]*GluTimeGap*1000;
opt.frameRate=5;
opt.colorbarLabel='∆F/F (%)';
opt.scaleBar=100/pixelsize;
opt.scaleBarText='100 µm';
opt.fontName='Arial';
opt.range=[0.01 0.1]*100;
writeMov4d(fullfile(fpath{f},'GTA_movieG_Syn#10'),mov2writeG_rgb,opt)

%% ===== 4. TRACE ANALYSIS — FIGURES ======================================
Syn2show = 7;               % synapse-of-interest to display (select here)

% ---- Select this synapse's precomputed peri-spike matrices ----
gspMat         = gspMat_all{Syn2show};
gspMat_voltage = gspMat_voltage_all{Syn2show};

% ---- Geodesic distances / ROI ordering relative to the selected synapse ----
interDendDist = zeros(1, nROI);
for i = 1:nROI
    interDendDist(i) = geodesic_distance(SkelDend, ...
        get_coord(AlignedVolt_ftprnt(:,:,i)), get_coord(GluResult.S_glu(:,:,Syn2show)));
end
[InterGludist_sorted, voltROIsorted] = sort(interDendDist, 'ascend');
GluFtprntDistFromSyn = zeros(1, size(GluResult.S_glu,3));
for si = 1:size(GluResult.S_glu,3)
    GluFtprntDistFromSyn(si) = geodesic_distance(SkelDend, Glu_coord(si,:), Glu_coord(Syn2show,:));
end
[~, GluROIsorted] = sort(GluFtprntDistFromSyn, 'ascend');

% ---- Classify this synapse's events as isolated / clustered ----
coincidence_threshold    = 1;
distancethreshold        = 50;
coactiveSyn_low_thresh   = 2;
coactiveSyn_upper_thresh = 5;
gspMat_coincidence = sum(gspMat(:,:,round(end/2)+[-coincidence_threshold:coincidence_threshold]),[3],'omitnan');
NcoactiveSyn   = sum(gspMat_coincidence(GluFtprntDistFromSyn*pixelsize<distancethreshold,:));
EventIsolated  = NcoactiveSyn <= coactiveSyn_low_thresh;
EventClustered = NcoactiveSyn >= coactiveSyn_upper_thresh;
IsolatedGTA    = squeeze(mean(gspMat_voltage(voltROIsorted,EventIsolated,:),2));
ClusteredGTA   = squeeze(mean(gspMat_voltage(voltROIsorted,EventClustered,:),2));

% ---- Figure 11: per-event glutamate raster + sorted voltage STA ----
gspMat2 = gspMat;
gspMat2(gspMat2==0) = NaN;
figure(11); clf; ax1 = [];
for ev = 1:size(gspMat,2)
    gspMat_sub        = squeeze(gspMat2(:,ev,:));
    gspMat_sub_sorted = squeeze(gspMat2(GluROIsorted,ev,:));
    [rr, cc] = find(~isnan(gspMat_sub_sorted));
    ax1 = [ax1 nexttile([1 1])];
    scatter(tax_glu_sub(cc)', rr, 20, 'k', '|', 'LineWidth', 2); hold all;
    scatter(tax_glu_sub, gspMat_sub(Syn2show,:), 'r|', 'LineWidth', 3);
    ylim([-1 50])
    yyaxis right;
    plot(tax_glu_sub, sum(~isnan(gspMat_sub),1), 'b-');
    ylim([0 15])
    ax1 = [ax1 nexttile([1 1])];
    Volt2show = squeeze(gspMat_voltage(voltROIsorted,ev,:));
    imagesc([-VoltLag:VoltLag]/1000, [1:nROI], movmean(Volt2show,15,2), caxis); colormap(turbo);
end
linkaxes(ax1,'x');
xlim([-150 150]/1000)

% ---- Figure 12: voltage decay for isolated vs clustered events ----
figure(12); clf; tiledlayout(2,4);

nexttile([1 1]);
histogram(NcoactiveSyn);
xlabel('# of coactive synapses'); ylabel('# of events'); box off;

nexttile([1 1]);
imagesc([-VoltLag:VoltLag]/1000, [1:nROI], IsolatedGTA, caxis); hold all; yyaxis right;
clusterActivity = squeeze(mean(gspMat(find(GluFtprntDistFromSyn*pixelsize<distancethreshold),find(EventIsolated),:),[1 2],'omitnan'));
plot(tax_glu_sub, clusterActivity, 'linewidth', 2); ylabel(sprintf('Fraction of coactive synapses within %d µm',distancethreshold));
title(sprintf(['Isolated glutamate-event \n triggered average voltage \n,' ...
    ' Coactive synapse ≤ %d within %d µm)'],coactiveSyn_low_thresh,distancethreshold));
xlim([-150 150]/1000); ylim([0 0.15]); xlabel('Time (sec)');

nexttile([1 1]);
imagesc([-VoltLag:VoltLag]/1000, [1:nROI], ClusteredGTA, caxis); hold all; yyaxis right;
title(sprintf(['Clustered glutamate-event \n triggered average voltage \n,' ...
    ' Coactive synapse ≥ %d within %d µm)'],coactiveSyn_upper_thresh,distancethreshold));
clusterActivity = squeeze(mean(gspMat(find(GluFtprntDistFromSyn*pixelsize<distancethreshold),find(EventClustered),:),[1 2],'omitnan'));
plot(tax_glu_sub, clusterActivity, 'linewidth', 2); ylim([0 0.15]);
ylabel(sprintf('Fraction of coactive synapses within %d µm',distancethreshold))
xlabel('Time (sec)'); xlim([-150 150]/1000);

nexttile([1 1]);
x_fit = [min(InterGludist_sorted):max(InterGludist_sorted)]'*pixelsize;
[y_fit_iso, t_consts_iso] = expfitDM_2(InterGludist_sorted'*pixelsize, mean(IsolatedGTA(:,VoltLag+[-30:30]),2),  x_fit, 100);
[y_fit_clu, t_consts_clu] = expfitDM_2(InterGludist_sorted'*pixelsize, mean(ClusteredGTA(:,VoltLag+[-30:30]),2), x_fit, 100);
plot(InterGludist_sorted*pixelsize, mean(IsolatedGTA(:,VoltLag+[-30:30]),2),  'color',[0 0.2 1]); hold all
plot(InterGludist_sorted*pixelsize, mean(ClusteredGTA(:,VoltLag+[-30:30]),2), 'color',[1 0.2 0]);
plot(x_fit, y_fit_iso, 'b--');
plot(x_fit, y_fit_clu, 'r--');
legend({sprintf('Isolated, Coactive synapse \n within %d µm ≤ %d, n = %d event',distancethreshold,coactiveSyn_low_thresh,sum(EventIsolated)),...
        sprintf('clustered, Coactive synapse \n within %d µm ≥ %d, n = %d event',distancethreshold,coactiveSyn_upper_thresh,sum(EventClustered)),...
        sprintf('Isolated fit, \\it{L} = %3.0f µm',t_consts_iso),...
        sprintf('Clusterd fit, \\it{L} = %3.0f µm',t_consts_clu)})
xlabel('Contour distance from synapse of interest (µm)');
ylabel('Voltage (-∆F/F)'); box off;

nexttile([1 4]);
show_footprnt_contour(AlignedVolt_ftprnt(:,:,voltROIsorted), VoltRefIm); hold all
scatter(Glu_coord(Syn2show,1), Glu_coord(Syn2show,2), 'ro')
drawScaleBar(100/pixelsize,'horizontal')

% ---- Figure 13: subthreshold depolarization vs cluster size ----
figure(13); Violin_wPoints(EventVoltageClustered, hsv(length(EventVoltageClustered)))
set(gca,'XTickLabel',[{'Isolated'}, arrayfun(@num2str, 2:9, 'UniformOutput', false)]);
xlabel('Cluster size (# synapse)'); ylabel('Subthreshold at glutamate event (∆F/F)');

% ---- Figure 14: isolated-event subthreshold, per synapse ----
figure(14); clf; tiledlayout(1,2);
ninetyprc_threshold = 5;
nexttile([1 1]);
scatter(IsolatedSubthreshold(:,2), IsolatedSubthreshold(:,3), 10, 'k', 'filled')
xlabel('90th percentile of cluster size'); ylabel(sprintf('Subthreshold at isolated \n glutamate event (∆F/F)'));
nexttile([1 1]);
IsolatedlikeSyn = IsolatedSubthreshold(:,2) < ninetyprc_threshold;
ClusterlikeSyn  = IsolatedSubthreshold(:,2) > ninetyprc_threshold;
p = Boxplot_wPoints2({IsolatedSubthreshold(IsolatedlikeSyn,3) IsolatedSubthreshold(ClusterlikeSyn,3)}, hsv(2));
drawPValueLines(p, 0, 'StepHeight', 0.005, 'TextYOffset', 0.002); box off;
ylim([-8 14]*0.001); ylabel(sprintf('Subthreshold at isolated \n glutamate event (∆F/F)'))
set(gca,'XTickLabel',{'Isolated-like synapse','Cluster-like synapse'})

%% ===== 5. MOVIE ANALYSIS — CALCULATIONS =================================


% ---- (c) Movie trigger times & isolated/clustered classification ----
Syn2look          = 7;
GluSpikeTime2look = find(GluSpike(Syn2look,:)>0);

N_pre  = 500;    % voltage frames BEFORE each trigger (ms at 1000 Hz)
M_post = 500;    % voltage frames AFTER  each trigger (ms at 1000 Hz)

output_STA_filename = fullfile(fpath{f}, 'GluTriggered_VoltSTA.mp4');
frame_rate = 30;
quality    = 95;

coincidence_threshold    = 1;
distancethreshold        = 50;
coactiveSyn_low_thresh   = 3;
coactiveSyn_upper_thresh = 6;

nWin         = N_pre + M_post + 1;      % total window length in volt frames
tau_axis     = (-N_pre : M_post);       % frame offset axis (0 = trigger)
tau_axis_glu = -ceil(N_pre/(exposuretime2*1000)) : ceil(M_post/(exposuretime2*1000));

glu_event_volt_frames = gluPeaks2VoltFrames(GluSpikeTime2look, GluResult.dFF_glu(Syn2look,:), ...
                                            GluResult.t_ax, VoltResult.t_ax, GluPeakWin);
[~, gspMat] = get_STA(GluSpike, GluSpike(Syn2show,:)>0, GluLag, GluLag);
[~, gspMat_voltage, validglu] = get_STA(VoltTrace_Norm, glu_event_volt_frames, N_pre, M_post);
validglu = ismember(glu_event_volt_frames, validglu);
GluSpikeTime2look     = GluSpikeTime2look(validglu);
glu_event_volt_frames = glu_event_volt_frames(validglu);
gspMat         = gspMat(:, validglu, :);
gspMat_voltage = gspMat_voltage - median(gspMat_voltage, 3);

gspMat_coincidence = sum(gspMat(:,:,round(end/2)+[-coincidence_threshold:coincidence_threshold]),[3],'omitnan');
NcoactiveSyn   = sum(gspMat_coincidence(GluFtprntDistFromSyn*pixelsize<distancethreshold,:));
EventIsolated  = NcoactiveSyn <= coactiveSyn_low_thresh;
EventClustered = NcoactiveSyn >= coactiveSyn_upper_thresh;
IsolatedGTA    = squeeze(mean(gspMat_voltage(voltROIsorted,EventIsolated,:),2));
ClusteredGTA   = squeeze(mean(gspMat_voltage(voltROIsorted,EventClustered,:),2));

% ---- (f) Normalize STA (all / isolated / clustered) ----
STA_volt               = zeros(sz_volt(1)*sz_volt(2), nWin);
STA_volt_sub           = zeros(sz_volt(1)*sz_volt(2), nWin);
STA_volt_sub_Isolated  = zeros(sz_volt(1)*sz_volt(2), nWin);
STA_volt_sub_Clustered = zeros(sz_volt(1)*sz_volt(2), nWin);

STA_volt(maskIdx,:)               = mean(TA_volt, 3, 'omitnan');
STA_volt_sub(maskIdx,:)           = mean(TA_volt_sub, 3, 'omitnan');
STA_volt_sub_Isolated(maskIdx,:)  = mean(TA_volt_sub(:,:,EventIsolated), 3, 'omitnan');
STA_volt_sub_Clustered(maskIdx,:) = mean(TA_volt_sub(:,:,EventClustered), 3, 'omitnan');

STA_volt               = toimg(STA_volt,               sz_volt(1), sz_volt(2));
STA_volt_sub           = toimg(STA_volt_sub,           sz_volt(1), sz_volt(2));
STA_volt_sub_Isolated  = toimg(STA_volt_sub_Isolated,  sz_volt(1), sz_volt(2));
STA_volt_sub_Clustered = toimg(STA_volt_sub_Clustered, sz_volt(1), sz_volt(2));

% ---- (g) Glutamate STA image, interpolated onto the voltage time grid ----
[~, TA_Glu] = get_STA(tovec(dFF_roi), GluSpikeTime2look, -tau_axis_glu(1), tau_axis_glu(end));
STA_Glu           = toimg(squeeze(mean(TA_Glu(:,:,:),2,'omitnan')),            sz_glu(1), sz_glu(2));
STA_Glu_Isolated  = toimg(squeeze(mean(TA_Glu(:,EventIsolated,:),2,'omitnan')), sz_glu(1), sz_glu(2));
STA_Glu_clustered = toimg(squeeze(mean(TA_Glu(:,EventClustered,:),2,'omitnan')),sz_glu(1), sz_glu(2));

TAxis_GluTrigger  = tau_axis_glu*exposuretime2;
TAxis_VoltTrigger = tau_axis*exposuretime1;
STA_Glu_interp           = interp1(TAxis_GluTrigger, tovec(STA_Glu)',           TAxis_VoltTrigger, 'linear', 'extrap').';
STA_Glu_Isolated_interp  = interp1(TAxis_GluTrigger, tovec(STA_Glu_Isolated)',  TAxis_VoltTrigger, 'linear', 'extrap').';
STA_Glu_clustered_interp = interp1(TAxis_GluTrigger, tovec(STA_Glu_clustered)', TAxis_VoltTrigger, 'linear', 'extrap').';

STA_Glu_interp           = toimg(STA_Glu_interp,           sz_glu(1), sz_glu(2));
STA_Glu_Isolated_interp  = toimg(STA_Glu_Isolated_interp,  sz_glu(1), sz_glu(2));
STA_Glu_clustered_interp = toimg(STA_Glu_clustered_interp, sz_glu(1), sz_glu(2));

STA_Glu_interp           = STA_Glu_interp           - median(STA_Glu_interp, 3);
STA_Glu_Isolated_interp  = STA_Glu_Isolated_interp  - median(STA_Glu_Isolated_interp(:,:,1:200), 3);
STA_Glu_clustered_interp = STA_Glu_clustered_interp - median(STA_Glu_clustered_interp(:,:,1:200), 3);

STA_Glu_interp(isnan(STA_Glu_interp))                     = 0;
STA_Glu_Isolated_interp(isnan(STA_Glu_Isolated_interp))   = 0;
STA_Glu_clustered_interp(isnan(STA_Glu_clustered_interp)) = 0;

% ---- (h) Polyline kymograph of the clustered STA movie ----
Vmeta = Device_Data{1,3};
Gmeta = Device_Data{1,4};
SynCoord2show      = Glu_coord(Syn2show,:);
SynCoord2show_volt = zeros(size(SynCoord2show));
Rvolt = refFromROI(size(VoltResult.ref_im),   double(Vmeta.ROI));
Rglu  = refFromROI(size(GluResult.AvgGluImg), double(Gmeta.ROI));
[xWorld, yWorld] = intrinsicToWorld(Rglu, SynCoord2show(:,1), SynCoord2show(:,2));
[SynCoord2show_volt(:,1), SynCoord2show_volt(:,2)] = transformPointsForward(tformReg, xWorld, yWorld);
[SynCoord2show_volt(:,1), SynCoord2show_volt(:,2)] = worldToIntrinsic(Rvolt, SynCoord2show_volt(:,1), SynCoord2show_volt(:,2));

[KymoTr, KymoROI] = polyLineKymo3(STA_volt_sub_Clustered, 10, 10, VoltResult.ref_im);
KymoROI_coord = cellfun(@(x) mean(x(1:end-1,:)), KymoROI, 'UniformOutput', false);
KymoROI_coord = cell2mat(KymoROI_coord');

Depol_glu  = mean(KymoTr(500+[-15:5],:), 1);
Depol_cmap = vec2cmap(Depol_glu, turbo(256));
ContourDist = [];
for roi = 1:length(KymoROI)
    ContourDist(roi) = geodesic_distance(DendriteMask_N1, KymoROI_coord(roi,:), SynCoord2show_volt);
end
[~, minROI] = min(ContourDist);

%% ===== 6. MOVIE ANALYSIS — FIGURES & MOVIE EXPORT =======================
% ---- Figure 22: kymograph summary ----
figure(22); clf; tiledlayout(2,2);
nexttile([1 2]);
imshow2(VoltResult.ref_im, []); hold all;
scatter(KymoROI_coord(:,1), KymoROI_coord(:,2), 25, Depol_cmap, 'filled');
l = scatter(SynCoord2show_volt(:,1), SynCoord2show_volt(:,2), 50, 'g^', 'filled'); hold all;
legend(l, 'Triggering synapse')

nexttile([1 1]);
imagesc([tau_axis], [1:length(KymoROI)], KymoTr'); hold all;
scatter([tau_axis(round(end/2))-145], minROI, 40, 'r>', 'filled');
xlim([-150 150]); xlabel('Time (ms)'); ylabel('ROI ID');

nexttile([1 1]);
plot(ContourDist(1:minROI)*pixelsize, Depol_glu(1:minROI)); hold all
plot(ContourDist(minROI:end)*pixelsize, Depol_glu(minROI:end));
xlabel('Contour distance (µm)'); ylabel('Subthreshold (∆F/F)'); box off;
legend({'Tip to synapse','Synapse to soma'})

% ---- Render 3x3 STA movie (subthreshold voltage | glutamate) ----
v = VideoWriter(output_STA_filename, 'MPEG-4');
v.FrameRate = frame_rate;
v.Quality   = quality;
open(v);
fig = figure('Color','k','Units','normalized','Position',[0.2 0.2 0.6 0.5]);
fprintf('Rendering STA movie...\n');
vcaxis = [-1 1]; gcaxis = [20 100];
OverlayScale = [1 0.5];
Str_rgb = grs2rgb(VoltResult.ref_im, gray(256), 20, 200);
bd = 5;
cmap_volt = gen_colormap([0 0.2 1; 0 0 0; 1 0 0], 256);
cmap_glu  = gen_colormap([0 0 0; 0 1 0], 256);

Alpha_thres     = 0.5;
Alpha_steepness = 10;
Alpha_kernal    = 4;

for fi = 1:2:length(tau_axis)

    alpha = 1 ./ (1 + exp(-Alpha_steepness * (abs(imgaussfilt(STA_volt_sub(bd:end-bd,bd:end-bd,fi),Alpha_kernal)) - Alpha_thres)));
    sub_rgb = grs2rgb(STA_volt_sub(bd:end-bd,bd:end-bd,fi),cmap_volt,vcaxis(1),vcaxis(2)).*alpha.*DendriteMask_N1(bd:end-bd,bd:end-bd);
    alpha2 = 1 ./ (1 + exp(-Alpha_steepness * (abs(imgaussfilt(STA_volt_sub_Isolated(bd:end-bd,bd:end-bd,fi),Alpha_kernal)) - Alpha_thres)));
    sub_iso_rgb = grs2rgb(STA_volt_sub_Isolated(bd:end-bd,bd:end-bd,fi),cmap_volt,vcaxis(1),vcaxis(2)).*alpha2.*DendriteMask_N1(bd:end-bd,bd:end-bd);
    alpha3 = 1 ./ (1 + exp(-Alpha_steepness * (abs(imgaussfilt(STA_volt_sub_Clustered(bd:end-bd,bd:end-bd,fi),Alpha_kernal)) - Alpha_thres)));
    sub_clu_rgb = grs2rgb(STA_volt_sub_Clustered(bd:end-bd,bd:end-bd,fi),cmap_volt,vcaxis(1),vcaxis(2)).*alpha3.*DendriteMask_N1(bd:end-bd,bd:end-bd);

    SGA_rgb     = transformCamera(Device_Data, tformReg, VoltResult.ref_im, STA_Glu_interp(:,:,fi));
    SGA_iso_rgb = transformCamera(Device_Data, tformReg, VoltResult.ref_im, STA_Glu_Isolated_interp(:,:,fi));
    SGA_clu_rgb = transformCamera(Device_Data, tformReg, VoltResult.ref_im, STA_Glu_clustered_interp(:,:,fi));

    SGA_rgb     = grs2rgb(SGA_rgb,     cmap_glu, gcaxis(1), gcaxis(2));
    SGA_iso_rgb = grs2rgb(SGA_iso_rgb, cmap_glu, gcaxis(1), gcaxis(2));
    SGA_clu_rgb = grs2rgb(SGA_clu_rgb, cmap_glu, gcaxis(1), gcaxis(2));

    clf(fig); set(fig, 'Color','k');
    tl = tiledlayout(fig, 3, 3, 'TileSpacing','compact', 'Padding','compact');

    ax7 = nexttile(tl, 1);
    imagesc(imgaussfilt(STA_volt_sub(bd:end-bd,bd:end-bd,fi),1), [vcaxis]); axis equal tight off;
    title(ax7, sprintf('Subthreshold voltage (%+.0f ms)', TAxis_VoltTrigger(fi)*1000), 'Color','w', 'FontSize',11);
    colormap(gray);

    ax7 = nexttile(tl, 4);
    imagesc(imgaussfilt(STA_volt_sub_Isolated(bd:end-bd,bd:end-bd,fi),1), [vcaxis]); axis equal tight off;
    colormap(gray);

    ax7 = nexttile(tl, 7);
    imagesc(imgaussfilt(STA_volt_sub_Clustered(bd:end-bd,bd:end-bd,fi),1), [vcaxis]); axis equal tight off;
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


