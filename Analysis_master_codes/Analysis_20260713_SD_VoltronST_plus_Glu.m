%% Combined Voltron-ST (soma voltage) + iGluSNFR (glutamate) extraction, spinning-disk data
%  Two-camera spinning-disk recording:
%    - Orange cam (Device_Data{3}, frames1 / mc*.bin , cam1_vsyn) -> Voltron-ST voltage
%    - Blue   cam (Device_Data{4}, frames2 / mc2*.bin, cam2_vsyn) -> iGluSNFR glutamate
%
%  Voltage soma: draw ROIs with clicky, find spikes on the rough ROI trace, build a
%  spike-triggered footprint pooled over all chunks, then project it back per chunk
%  with extractVoltronST (project mode). Streams one chunk at a time (low memory).
%  Glutamate is extracted with extractGluSNFR3 (per-tile detect -> merge -> project),
%  same as Analysis_20260206_iGluSNFRextraction.m.
%
%  Assumes motion-corrected bins already exist:
%    voltage : mc%02d.bin      + mcTrace%02d.mat
%    glut.   : mc2%02d.bin      + mc2Trace%02d.mat
%
%  Outputs per file:  Volt_Result.mat , Glu_Result.mat

%% Section 0 : Setup
clear
clc;

% Remove any shadow of the built-in 'graph' (old cvx/sdpt3 defines graph.m,
% which breaks the centroid clustering in extractVoltronST / extractGluSNFR3).
gshadow = which('graph','-all');
for gi = 1:numel(gshadow)
    if ~contains(gshadow{gi}, fullfile('toolbox','matlab'))   % keep MATLAB's built-in
        rmpath(fileparts(gshadow{gi}));
        fprintf('Removed shadowing graph.m folder from path: %s\n', fileparts(gshadow{gi}));
    end
end

[~, ~, raw] = xlsread('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement.xlsx', 'Sheet1', 'C5:AA31');
load('/Volumes/BHL18TB_D2/20260203_SD_V2+iGluSNFR4/25X_transformationMatrix.mat');   % tformReg
fpath           = raw(:,1)';
V2moviemaxTime  = 15000;   % voltage movie chunk length (frames)
GlumoviemaxTime = 5000;    % glutamate movie chunk length (frames)
o_laser         = 200001;  % 607/488 modulation onset (DAQ sample index, ~2 s after acq start)
mod488          = 200001;
foi             = [14:19];   % files of interest
set(0,'DefaultFigureWindowStyle','docked')

%% Section 1a : Voltage — timing setup + reference-chunk ROI (clicky)
%  Manual workflow: draw soma ROIs on a reference chunk (clicky). Section 1b then
%  builds the rough soma trace across all chunks and finds spikes; Section 1c builds
%  a spike-triggered footprint (pooled over all chunks) and projects it back.
%  Per-chunk residual movies come from local function loadVoltChunk (end of file).
for f=[14:19]
refSeg = 1;     % chunk to draw ROIs on
fprintf('=== Voltage (Voltron-ST) file #%2.0f : %s\n', f, fpath{f});
Result = struct();
load(fullfile(fpath{f}, 'output_data.mat'));   % Device_Data

sz       = double(Device_Data{3}.ROI([2 4]));     % [nCol nRow] orange cam
DAQ_rate = Device_Data{1,2}.Counter_Inputs(1,1).rate;

% ---- Reconstruct voltage camera clock (fill tail with synthetic triggers) ----
cam1_vsyn = Device_Data{1,2}.Counter_Inputs(1,1).data;
start_idx = min(find(cam1_vsyn == max(cam1_vsyn)));
end_idx   = length(Device_Data{1,2}.buffered_tasks(1,2).channels(1,1).data);
CamTrigger1_DAQax = find((cam1_vsyn(2:end) - cam1_vsyn(1:end-1)) > 0);
segment_size = round(median(diff(CamTrigger1_DAQax)));
last_val   = cam1_vsyn(start_idx - 1);
n_to_add   = end_idx - start_idx + 1;
n_segments = ceil(n_to_add / segment_size);
added_part = repelem((last_val + 1 : last_val + n_segments)', segment_size);
added_part = added_part(1:n_to_add);
cam1_vsyn(start_idx:end_idx) = added_part;

% ---- Movie segmentation + voltage time axis ----
VmovTimesegments = (cam1_vsyn(o_laser)+2) : V2moviemaxTime : cam1_vsyn(end);
VmovTimesegments(end+1) = cam1_vsyn(end);
nFrame2analyze = VmovTimesegments(end) - VmovTimesegments(1);
CamTrigger1_DAQaxVec = ind2vec(length(cam1_vsyn), CamTrigger1_DAQax(1):segment_size:cam1_vsyn(end)*segment_size, 1);
CamTrigger1_DAQaxVec = CamTrigger1_DAQaxVec(o_laser:end);
FirstFrameDAQax = find(CamTrigger1_DAQaxVec > 0, 1);
t_vol = (FirstFrameDAQax : segment_size : (FirstFrameDAQax + segment_size*nFrame2analyze)) / DAQ_rate;

nSeg        = length(VmovTimesegments)-1;
readFrames  = diff(VmovTimesegments);   % frames requested per chunk

Result.VmovTimesegments=VmovTimesegments;
Result.t_ax = t_vol;

% ---- Draw ROIs on the reference chunk ----
[movR_ref, mov_mc_ref] = loadVoltChunk(fpath{f}, refSeg, sz, readFrames(refSeg));
Result.ref_im   = mean(mov_mc_ref, 3);
[roi_points, ~] = clicky(movR_ref, Result.ref_im);   % click somata, right-click to finish
Result.ROIpoly  = roi_points;
nROI = numel(roi_points);
fprintf('  %d ROI(s) drawn on chunk %d\n', nROI, refSeg);
clear movR_ref mov_mc_ref
save(fullfile(fpath{f},'VoltResult.mat'),'Result','-v7.3');
fprintf('File #%d, ROI saved\n',f);
end
 %% Section 1b : Voltage — rough soma trace (all chunks) + spike finding
%  Re-apply the drawn ROIs to every chunk to build a full-length rough trace,
%  then detect spikes (used only to build the footprint in 1c).

staTau   = -20:20;    % frames averaged around each spike to form the footprint
dilatePx = 5;       % px dilation of each ROI bounding its footprint
[H,W]    = size(Result.ref_im);

for f=[15:19]
load(fullfile(fpath{f},'VoltResult.mat'));
load(fullfile(fpath{f}, 'output_data.mat'));   % Device_Data

nROI = numel(Result.ROIpoly);
sz       = double(Device_Data{3}.ROI([2 4]));     % [nCol nRow] orange cam
readFrames  = diff(Result.VmovTimesegments);
nSeg        = length(Result.VmovTimesegments)-1;
rough_tr  = [];
segFrames = zeros(nSeg,1);
Result.mc = [];
for j = 1:nSeg
    [movR, ~, mtr, nT] = loadVoltChunk(fpath{f}, j, sz, readFrames(j));
    segFrames(j) = nT;
    Result.mc    = [Result.mc; mtr];
    intens   = apply_clicky(Result.ROIpoly, movR, 'no');   % nT x nROI (spike-positive)
    rough_tr = [rough_tr, intens'];                    % nROI x Ttotal
    fprintf('  rough trace chunk %2.0d/%2.0d\n', j, nSeg);
end

% z-score (high-pass) and detect spikes per ROI
tr_hi = rough_tr - movmedian(rough_tr, 50, 2);
tr_z  = tr_hi ./ get_threshold(tr_hi, 1);
sp    = find_spike_bh(tr_z, 2, 1);                      % nROI x Ttotal binary
sp(:,round(end*0.8):end)=0;
Result.spike_rough = sp;
Result.rough_tr = rough_tr;
fprintf('  spikes per ROI: %s\n', mat2str(sum(sp,2)'));

% quick look
figure(3*f-1); clf;
show_traces_spikes(rescale2(rough_tr(:,1:580000),2),sp(:,1:580000),sum(sp(:,1:580000),1));
drawnow;

% Section 1c : Voltage — spike-triggered footprint (pooled) + projection + save

% ---- PASS A: accumulate spike-triggered footprint over ALL chunks (streaming) ----
% get_STA returns the per-chunk MEAN window (pixels x nTau); scale by the chunk
% spike count and divide by the pooled count later to get the global mean STA.
tau1 = -staTau(1); tau2 = staTau(end);   % window = [-tau1:tau2]
staSum = zeros(H,W,length(staTau),nROI);
staCnt = zeros(1,nROI);
for j = 1:nSeg
    off = sum(segFrames(1:j-1));                         % global frame offset of this chunk
    if ~any(sp(:, off+1 : off+segFrames(j)), 'all')      % no spikes in this chunk -> skip loading
        fprintf('  STA chunk %2.0d/%2.0d skipped (no spikes)\n', j, nSeg);
        continue
    end
    movR   = loadVoltChunk(fpath{f}, j, sz, readFrames(j));
    movVec = tovec(movR);
    for n = 1:nROI
        [~, mat2cat, tp] = get_STA(movVec, sp(n, off+1 : off+segFrames(j)), tau1, tau2);
        if isempty(tp), continue; end                   % no spikes here (edge spikes dropped inside get_STA)
        staSum(:,:,:,n) = staSum(:,:,:,n) + toimg(permute(sum(mat2cat,2),[1 3 2]), H, W);
        staCnt(n)       = staCnt(n) + numel(tp);
    end
    fprintf('  STA chunk %2.0d/%2.0d\n', j, nSeg);
end

% ---- Build S_full: mask each footprint to its ROI vicinity, then L2-normalize ----
S_full = zeros(H,W,nROI,'double');
for n = 1:nROI
    if staCnt(n) > 0
        fp = staSum(:,:,:,n) / staCnt(n);
        fp = fp - median(fp(:,:,1:15),3);
        fp = imgaussfilt(max(fp ,[],3),1);
    else
        warning('ROI %d has no spikes; footprint left empty.', n);
        fp = zeros(H,W);
    end
    m  = imdilate(poly2mask(Result.ROIpoly{n}(:,1), Result.ROIpoly{n}(:,2), H, W), strel('disk',dilatePx));
    fp = fp .* m;
    S_full(:,:,n) = fp / (sqrt(sum(fp(:).^2)) + eps);
end
Result.ftprnt = S_full;
Result.STAmovieRough=staSum;
Result.STAcountRough=staCnt;

% ---- PASS B: project the footprint back onto each chunk (streaming) ----
optsP = struct(); optsP.mode = 'project';
T_all = []; F0_all = [];
for j = 1:nSeg
    [movR, mov_mc] = loadVoltChunk(fpath{f}, j, sz, readFrames(j));
    tmp = extractVoltronST(movR, mov_mc, optsP, S_full);   % movR is spike-positive
    T_all  = [T_all,  tmp.T_glu];
    F0_all = [F0_all, tmp.F0_glu];
    fprintf('  project chunk %2.0d/%2.0d\n', j, nSeg);
end

% ---- Assemble result (global dF/F baseline across the whole recording) ----

Result.F0      = F0_all;
Result.traces           = T_all;

% ---- Diagnostics + save ----
figure(2*f); clf; tiledlayout(2,2);
ax1 = nexttile([1 1]); show_footprnt_contour(Result.ftprnt, Result.ref_im); colormap(ax1,gray(256));
nexttile([1 1]); imshow2(im_merge(Result.ftprnt, lines(size(Result.ftprnt,3)))*10, []);
nexttile(3,[1 2]); plot(rescale2(Result.traces,2)' + (1:size(Result.traces,1)));
sgtitle(fpath{f}, 'Interpreter','none'); drawnow;
saveas(gcf, fullfile(fpath{f}, 'Volt_extraction_result.fig'));
saveas(gcf, fullfile(fpath{f}, 'Volt_extraction_result.png'), 'png');
save(fullfile(fpath{f}, 'Volt_Result.mat'), 'Result', '-v7.3');
fprintf('Saved Volt_Result.mat to %s\n', fpath{f});

end
%% Section 2 : Glutamate extraction via extractGluSNFR3
%  Interactive per file (PC-regression prompt + optional UIs inside extractGluSNFR3).
for f=[14:19]
load(fullfile(fpath{f}, 'output_data.mat'));   % Device_Data
DAQ_rate = Device_Data{1,2}.Counter_Inputs(1,1).rate;

% orange cam (voltage) — for co-registration image only
nCol1 = double(Device_Data{3}.ROI(2));
nRow1 = double(Device_Data{3}.ROI(4));
exposuretime1 = Device_Data{3}.exposuretime;

% blue cam (glutamate)
nCol2 = double(Device_Data{4}.ROI(2));
nRow2 = double(Device_Data{4}.ROI(4));
exposuretime2 = Device_Data{4}.exposuretime;

% ---- Reconstruct glutamate (blue) camera clock ----
cam2_vsyn = Device_Data{1,2}.Counter_Inputs(1,2).data;
cam2_vsyn_trig = find(diff(cam2_vsyn) == 1) + 1;
segment_size2  = cam2_vsyn_trig(10) - cam2_vsyn_trig(9);
start_idx = min(find(cam2_vsyn == max(cam2_vsyn)));
end_idx   = length(Device_Data{1,2}.buffered_tasks(1,2).channels(1,2).data);
last_val   = cam2_vsyn(start_idx - 1);
n_to_add   = end_idx - start_idx + 1;
n_segments = ceil(n_to_add / segment_size2);
added_part = repelem((last_val + 1 : last_val + n_segments)', segment_size2);
added_part = added_part(1:n_to_add);
cam2_vsyn(start_idx:end_idx) = added_part;

cam2_trig = find(diff(cam2_vsyn) == 1) + 1;                       % DAQ frame index
t_start   = cam2_trig(cam2_vsyn(mod488)+2) - mod488 + 1;
if cam2_vsyn(end) < GlumoviemaxTime
    GluemovTimesegments = [(cam2_vsyn(mod488)+2) cam2_vsyn(end)];
else
    GluemovTimesegments = (cam2_vsyn(mod488)+2) : GlumoviemaxTime : cam2_vsyn(end);
    GluemovTimesegments(end+1) = cam2_vsyn(end);
end
nFrame2analyze2 = GluemovTimesegments(end) - GluemovTimesegments(1);
t_glu = (t_start : segment_size2 : (t_start + segment_size2*nFrame2analyze2)) / DAQ_rate;

% voltage average image for co-registration
Vmov_read_tmp = double(readBinMov_times([fpath{f} '/mc' num2str(1,'%02d') '.bin'], nRow1, nCol1, [2000:3000]));
AvgVoltageImg = mean(Vmov_read_tmp, 3);

Ntile = length(GluemovTimesegments) - 1;
T     = GluemovTimesegments(end) - GluemovTimesegments(1);
GluResult = struct();

% ---------- PASS 1: per-tile detection ----------
tileResults        = cell(Ntile,1);
RegressComponentTile = cell(Ntile,1);
mcTrace_glu = [];
opts = struct(); opts.doPlot = true;

for k = 1:Ntile
    fprintf('Glu detect: reading tile %1.0f / %1.0f ...\n', k, Ntile);
    mov2_mc = double(readBinMov_times([fpath{f} '/mc2' num2str(k,'%02d') '.bin'], nRow2, nCol2, ...
        [1:(GluemovTimesegments(k+1)-GluemovTimesegments(k))]));
    load([fpath{f} '/mc2Trace' num2str(k,'%02d') '.mat']);

    meanF = squeeze(mean(mov2_mc,[1 2]));
    y_fit = expfitDM_2((1:size(mov2_mc,3))', meanF, (1:size(mov2_mc,3))', [10000 100]);
    mov2_mc_filt    = mov2_mc ./ reshape(y_fit,1,1,[]);
    mov2_mc_filt_lw = movmedian(mov2_mc_filt, round(5/exposuretime2), 3);
    mov_sub = mov2_mc_filt - mov2_mc_filt_lw;

    [V, D, u_sub] = get_eigvector(tovec(imgaussfilt(mov_sub,1))', 40);
    u_sub = reshape(u_sub, nRow2(1), nCol2(1), []);

    if ~isfield(GluResult,'bvMask')
        [~, GluResult.bvMask] = get_ROI(mean(mov2_mc,3));
    end
    figure(2); clf; imshow2_patch(u_sub); drawnow;
    n = input("PCs to regress out\n");
    RegressComponentTile{k} = n;
    GluResult.RegressComponentTile{k} = n;

    mov_sub = SeeResiduals(mov_sub, V(:,n));
    mov_sub = SeeResiduals(mov_sub, mcTrace2.xymean(1:(GluemovTimesegments(k+1)-GluemovTimesegments(k)),:));
    mcTrace_glu = [mcTrace_glu; mcTrace2.xymean(1:(GluemovTimesegments(k+1)-GluemovTimesegments(k)),:)];
    mov2_mc_filt2 = mov_sub + mov2_mc_filt_lw;
    mov_sub = mov_sub .* (max(GluResult.bvMask,[],3) == 0);

    if ~isfield(tileResults{1},'opts')   % first tile: set detection opts
        opts.diagSaveDir     = fpath{f};
        opts.diagPrefix      = 'exp01';
        opts.diagZoom_pad    = 15;
        opts.diagMaxSyn      = 16;
        opts.eventMaxArea    = 30;
        opts.nmfReps         = 5;
        opts.nmfDiagMaxEvent = 6;
        opts.showNMFDiag     = true;
    else
        opts = tileResults{1}.opts;
        opts.uiClusterSize = 0;
        opts.plotCentroids = 0;
        opts.uiSynFilter   = 1;
        opts.uiSeedThresh  = 1;
    end

    tileResults{k} = extractGluSNFR3(mov_sub, mov2_mc_filt2, exposuretime2, opts);
    drawnow;
end

% ---------- Merge synapse dictionary across tiles ----------
if Ntile > 1
    Dict   = mergeSynapseDict(tileResults, 3);
    S_full = Dict.S_glu;

    opts2 = struct(); opts2.mode = 'project';
    Rproj = cell(Ntile,1);
    for k = 1:Ntile
        fprintf('Glu project: reading tile %1.0f / %1.0f ...\n', k, Ntile);
        try
            mov2_mc = double(readBinMov([fpath{f} '/mc2' num2str(k,'%02d') '.bin'], nRow2, nCol2));
        catch
            mov2_mc = double(readBinMov_times([fpath{f} '/mc2' num2str(k,'%02d') '.bin'], nRow2, nCol2, ...
                [1:(GluemovTimesegments(k+1)-GluemovTimesegments(k))]));
        end
        meanF = squeeze(mean(mov2_mc,[1 2]));
        y_fit = expfitDM_2((1:size(mov2_mc,3))', meanF, (1:size(mov2_mc,3))', [10000 100]);
        mov2_mc_filt    = mov2_mc ./ reshape(y_fit,1,1,[]);
        mov2_mc_filt_lw = movmedian(mov2_mc_filt, 100, 3);
        mov_sub = mov2_mc_filt - mov2_mc_filt_lw;
        [V, D, u_sub] = get_eigvector(tovec(imgaussfilt(mov_sub,1))', 40);
        mov_sub = SeeResiduals(mov_sub, V(:,RegressComponentTile{k}));
        mov2_mc_filt2 = mov_sub + mov2_mc_filt_lw;
        Rproj{k} = extractGluSNFR3(mov_sub, mov2_mc_filt2, exposuretime2, opts2, S_full);
    end
else
    Rproj  = tileResults;
    S_full = Rproj{1,1}.S_glu;
end

frameRate = 1/exposuretime2;
Nsyn = size(S_full,3);

% ---------- Concatenate tiles along time ----------
T_all = []; F0_all = []; dFF_all = [];
for k = 1:numel(Rproj)
    R = Rproj{k};
    if isempty(R) || isempty(R.T_glu), continue; end
    assert(size(R.T_glu,1) == Nsyn, 'Tile %d has wrong Nsyn.', k);
    T_all = [T_all, R.T_glu];
    if isfield(R,'F0_glu')  && ~isempty(R.F0_glu),  F0_all  = [F0_all,  R.F0_glu];  end
    if isfield(R,'dFF_glu') && ~isempty(R.dFF_glu), dFF_all = [dFF_all, R.dFF_glu]; end
end

GluResult.mc        = mcTrace_glu;
GluResult.S_glu     = S_full;
GluResult.T_glu     = T_all;
GluResult.F0_glu    = F0_all;
GluResult.dFF_glu   = dFF_all;
GluResult.frameRate = frameRate;
GluResult.AvgGluImg = mean(mov2_mc,3);
GluResult.t_ax      = t_glu;

save(fullfile(fpath{f}, 'Glu_Result.mat'), '-struct', 'GluResult', '-v7.3');
fprintf('Saved Glu_Result.mat to %s\n', fpath{f});
end
%% Section 3 : Combined co-registration + overlay (voltage soma + glutamate synapses)
%  Requires Volt_Result.mat (Section 1) and Glu_Result.mat (Section 2) for file f.
VoltResult = importdata(fullfile(fpath{f}, 'Volt_Result.mat'));
GluResult  = importdata(fullfile(fpath{f}, 'Glu_Result.mat'));

S_full  = GluResult.S_glu;
dFF_all = GluResult.dFF_glu;
frameRate = GluResult.frameRate;

% Map glutamate footprints into the voltage camera frame (and vice versa)
[alignedGlu, Rglu, Rvolt] = transformCamera(Device_Data, tformReg, VoltResult.ref_im, max(S_full,[],3));
alignedVolt = transformCamera_O2B(Device_Data, tformReg, VoltResult.ref_im, GluResult.AvgGluImg);

% Keep only synapses that fall on the voltage-labelled soma/dendrite
glucoord = round(get_coord(S_full));
validSynapse = zeros(1, size(glucoord,1));
for s = 1:size(glucoord,1)
    validSynapse(s) = alignedVolt(glucoord(s,2), glucoord(s,1)) > 0;
end

% dF/F traces of valid synapses
figure(99); clf;
line_color = hsv(sum(validSynapse));
l = plot(GluResult.t_ax, dFF_all(validSynapse>0,:)' + (1:sum(validSynapse))*1);
arrayfun(@(l,c) set(l,'Color',c{:}), l, num2cell(line_color,2));
title('Glutamate dF/F (synapses on voltage cell)');
xlabel('Time (s)'); ylabel('Spine ID'); axis tight;
saveas(gca, fullfile(fpath{f}, 'Glu_trace.fig'));

% Footprint overlays
figure(100); clf; tiledlayout(1,3);
offsetscale = 0.7;
coloredS_glu = showScaleImage(S_full(:,:,validSynapse>0), 1:sum(validSynapse), 'hsv');

nexttile([1 1]);
imshow2(coloredS_glu + grs2rgb(alignedVolt, gray(256))*offsetscale, []);
title('Glu synapses (color) + voltage (gray)');

nexttile([1 1]);
imshow2(im_merge(cat(3, mat2gray(alignedVolt), mat2gray(max(S_full,[],3))), [0 0.6 1; 1 0 0]), []);
title('Voltage (blue) + Glu (red), glu frame');

nexttile([1 1]);
imshow2(im_merge(cat(3, mat2gray(GluResult.AvgGluImg), mat2gray(max(S_full,[],3))), [0 0.6 1; 1 0 0]), []);
title('Glu mean (blue) + Glu footprints (red)');
saveas(gca, fullfile(fpath{f}, 'Glu_Volt_footprint.png'));