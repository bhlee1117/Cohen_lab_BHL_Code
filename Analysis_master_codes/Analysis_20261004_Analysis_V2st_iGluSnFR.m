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
iswindow=0;
try
    [~, ~, raw] = xlsread('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement_V2.xlsx', 'Sheet1', 'C5:AA45');
catch
    raw = readcell(macToWindowsPath('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Prism_V2+Glu_Data_Arrangement_V2.xlsx'),...
        'Sheet', 'Sheet1', 'Range', 'C5:AA45');
    iswindow=1;
end

fpath           = raw(:,1)';
VRfpath           = raw(:,8)';
VRbw= (raw(:,9));
fpath_valid=cell2mat(cellfun(@(x) any(~ismissing(x)),fpath,'UniformOutput',false));
VRfpath_valid=cell2mat(cellfun(@(x) any(~ismissing(x)),VRfpath,'UniformOutput',false));
if iswindow
    for f=find(fpath_valid>0)
        fpath{f}=macToWindowsPath(fpath{f});        
    end
    for f=find(VRfpath_valid>0)
        VRfpath{f}=macToWindowsPath(VRfpath{f});        
    end
end
V2moviemaxTime  = 15000;   % voltage movie chunk length (frames)
GlumoviemaxTime = 3500;    % glutamate movie chunk length (frames)
o_laser         = 200001;  % 607/488 modulation onset (DAQ sample index, ~2 s after acq start)
mod488          = 200001;
foi             = [15 16 18];   % files of interest
set(0,'DefaultFigureWindowStyle','docked')

% --- extra columns + transform used by the f=3/f=5 sections below (cross-platform) ---
StructureData = raw(:,6);            % NB: paths in here also need macToWindowsPath on Windows
tformFile = "/Users/bhlee1117/Documents/GitHub/Cohen_lab_BHL_Code/16X_transformationMatrix.mat";
if iswindow, tformFile = macToWindowsPath(tformFile); end
load(tformFile);                     % tformReg

%% Section 1 : V2-ST / iGluSnFR expression overlap in the structure z-stack
%  Uses Fullset_Result.mat from Analysis_20261002_SD_VoltronST_plus_Glu.m (Section 5). The structure
%  stack is cut to the recording FOV (glutamate camera pixels: Structure.rows/cols) and analysed in
%  the overlapping field: the part of the glutamate FOV the voltage camera also sees
%  (Volt.ref_im_glu > 0) and the structure covers. Both channels are the max over the matched
%  slice +-ov.zHalf.
%   Orange (V2-ST)  : each voltage ROI's footprint on the glutamate image (Volt.ftprnt_glu) marks
%                     where the neuron is; its soma is the brightest soma-sized orange blob inside
%                     the footprint. When the stack shows none (soma out of the slab or too dim),
%                     it is searched the same way in the voltage reference image (Volt.ref_im_glu),
%                     brought to the structure's orange brightness first (background and the
%                     orange somata found in the stack matched); that neuron's orange is then
%                     measured in this image. Still none: a disk at the footprint peak.
%   Blue (iGluSnFR) : somata from a soma-scale Hessian blob detector + seeded watershed.
%   A blue soma that overlaps an orange soma is the same neuron and keeps the blue boundary.
%   Every neuron's ROI is measured in both channels against a ring around it:
%     z = (mean in ROI - median of ring) / robust SD of ring;  z >= ov.posZ = expressing.
%  Before that, Volt.ref_im_glu (reached the glutamate frame through tformReg) is registered to the
%  stack's orange channel (reached it through the glutamate mean image) by a translation, and the
%  footprints are shifted with it.
%  A neuron seen in several recordings of the same structure stack is counted once (dup flag).
ov.foi         = 14:29;       % recordings with Fullset_Result.mat (those without an orange stack are skipped)
ov.zHalf       = 1;           % slices on each side of the matched slice (max projection)
ov.edge_px     = 3;           % erosion of the voltage-camera FOV edge
ov.hpSigma     = 30;          % background (Gaussian sigma, px) removed before soma detection
ov.somaSigma   = 4;           % soma scale of the blob detector (px; soma 12-15 um = 14-17 px)
ov.seedZ       = 4;           % seed: smoothed high-passed peak above this many robust SD of the field
ov.seedMinDist = 8;           % a seed closer than this to a stronger one is dropped (px)
ov.edgeFrac    = 0.3;         % soma edge: local background + the lower of edgeFrac * (peak - background)
ov.edgeZ       = 3;           %   and edgeZ * background SD (a bright soma is not cut down to its core);
ov.edgeFracMax = 0.6;         %   raised up to edgeFracMax while the soma is too large or not compact
ov.wsSigma     = 2.5;         % smoothing of the watershed image (px), bridges the dark nucleus
ov.mergeFrac   = 0.3;         % two seeds are one soma when the profile between them stays above
                              %   mergeFrac of the weaker peak (iGluSnFR somata have a dark nucleus)
ov.mergeRatio  = 0.3;         %   and the weaker peak is at least mergeRatio of the stronger
ov.openR       = 3;           % opening radius that strips dendrites from a soma (px)
ov.maxR        = 18;          % soma pixels at most this far from its seed(s) (px)
ov.area        = [100 1000];  % soma area range (px; 100 px = 11 um diameter: rejects puncta and fragments)
ov.edgeCut     = 0.05;        % a soma whose 2-px rim outside the field exceeds this fraction of its area is dropped
ov.minSolidity = 0.6;         % compactness (pyramidal somata with processes are ~0.65)
ov.ftThr       = 0.5;         % footprint region: smoothed footprint > ftThr * its max
ov.refMaxShift = 30;          % fine registration of Volt.ref_im_glu onto the stack's orange: search +-px
ov.refMinCorr  = 0.3;         %   shift applied when the correlation peak reaches this
ov.refSeedZ    = 1.5;         % seed threshold in the voltage reference image: the search is confined to the
                              %   footprint of a drawn ROI, and the mean image has low soma contrast
ov.calMinPix   = 60;          % ref image brightness matched on the stack's orange somata when they
                              %   cover this many px, else on the 99.5th percentiles
ov.somaR       = 6;           % disk radius when no orange soma is found in the footprint (px)
ov.matchFrac   = 0.3;         % same neuron: overlap / smaller ROI area
ov.bgRing      = [3 12];      % background ring, px outside the ROI
ov.posZ        = 2;           % expressing when z >= posZ (each channel)
ov.dedupDist   = 8;           % same neuron in another recording of the same stack: within this
ov.dedupDz     = 3;           %   many px and slices
ov.maxAspect   = Inf;         % soma long/short axis limit (Section 1: footprints already locate the soma)
ov.wsMaxAspect = 2.5;         % Section 2 (no voltage ROIs): rejects dendrite segments
ov.wsSeedZ     = [4 10];      % Section 2 seed threshold [blue orange]; the background-subtracted orange
                              %   stack has little noise, so faint blobs in dark areas pass at 4
ov.pixel_um    = 0.887;       % structure camera um/px (scale bar)
ov.saveDir     = '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_V2+Glu/';
if iswindow, ov.saveDir = macToWindowsPath(ov.saveDir); end

Ov = {};
for f = ov.foi
    ff = fullfile(fpath{f}, 'Fullset_Result.mat');
    if ~isfile(ff), fprintf('#%d no Fullset_Result.mat, skipped\n', f); continue; end
    S = load(ff, 'Fullset'); F = S.Fullset;
    if ~isfield(F.Structure, 'cropOrange'), fprintf('#%d structure stack has no orange channel, skipped\n', f); continue; end
    R = expressionOverlap(F.Volt, F.Structure, size(F.Glu.AvgGluImg), ov);
    R.file = f; R.fpath = fpath{f}; R.mouse = F.mouse; R.structure = F.Structure.name;
    Ov{end+1} = R; %#ok<SAGROW>
    fprintf('#%d %s: %d neurons in the overlapping field (%d blue somata, %d voltage ROIs, %d matched)\n', f, ...
        R.structure, numel(R.neurons), numel(R.blue), numel(R.orange), nnz([R.neurons.voltROI] > 0 & [R.neurons.blueCell] > 0));
    fprintf('    ref_im_glu -> stack orange: shift [%+d %+d] px, corr %.2f (unshifted %.2f)%s\n', R.refReg.dx, R.refReg.dy, ...
        R.refReg.corr, R.refReg.corr0, repmat(' NOT APPLIED', 1, ~R.refReg.applied));
end
Ov = [Ov{:}];
T = overlapTable(Ov, ov);
U = T(~T.dup, :);                                        % unique neurons
printOverlapStats(U, T);
if ~isfolder(ov.saveDir), mkdir(ov.saveDir); end
save(fullfile(ov.saveDir, 'ExpressionOverlap_V2_iGluSnFR.mat'), 'Ov', 'T', 'ov', '-v7.3');

%% Section 1b : Representative structure image (overlapping field)
ov.rep = 16;                                             % representative recording
R = Ov([Ov.file] == ov.rep);
fig = overlapRepFigure(R, ov, 4000+ov.rep);
savePNG(fig, fullfile(ov.saveDir, sprintf('ExpressionOverlap_rep_%02d.png', ov.rep)));

%% Section 1c : Expression overlap statistics
fig = overlapStatsFigure(T, ov, 4100);
savePNG(fig, fullfile(ov.saveDir, 'ExpressionOverlap_stats.png'));

%% Section 2 : Expression overlap over the whole structure image
%  As Section 1, but over the whole stack slab of each recording (matched slice +-ov.zHalf of the
%  full z-stack, not only the recording FOV). There are no voltage ROIs to guide the orange search
%  here, so orange somata are detected in the stack the same way as blue somata (with an elongation
%  limit and a higher orange seed threshold, ov.wsMaxAspect / ov.wsSeedZ), and an orange and a blue
%  soma that overlap are one neuron (blue boundary). A neuron within ov.dedupDist of a
%  Section 1 voltage-recorded neuron gets its voltage ROI number (recorded); inFOV marks neurons
%  inside the recording FOV. Recordings imaging the same stack slab are counted once (dup).
Ws = {};
for k = 1:numel(Ov)
    S = load(fullfile(Ov(k).fpath, 'Structure_Result.mat'), 'Structure'); St = S.Structure;
    if iswindow, St.blueFile = macToWindowsPath(St.blueFile); St.orangeFile = macToWindowsPath(St.orangeFile); end
    W2 = wholeStackOverlap(St, Ov(k), ov);
    W2.file = Ov(k).file; W2.fpath = Ov(k).fpath; W2.mouse = Ov(k).mouse; W2.structure = Ov(k).structure;
    Ws{end+1} = W2; %#ok<SAGROW>
    fprintf('#%d %s: %d neurons in the stack slab (%d blue, %d orange somata, %d both; %d in the recording FOV, %d recorded)\n', ...
        W2.file, W2.structure, numel(W2.neurons), numel(W2.blue), numel(W2.orangeCells), ...
        nnz(strcmp({W2.neurons.source}, 'blue+orange')), nnz([W2.neurons.inFOV]), nnz([W2.neurons.voltROI] > 0));
end
Ws = [Ws{:}];
T2 = overlapTable(Ws, ov);
U2 = T2(~T2.dup, :);
printOverlapStats(U2, T2, 'whole structure image');
fprintf('  inside the recording FOV: %d neurons, outside: %d\n', nnz(U2.inFOV), nnz(~U2.inFOV));
WsSave = rmfield(Ws, {'B', 'O', 'valid'});               % stack images are re-read for figures
save(fullfile(ov.saveDir, 'ExpressionOverlap_V2_iGluSnFR_wholeStack.mat'), 'WsSave', 'T2', 'ov', '-v7.3');

%% Section 2b : Representative whole structure image
W2 = Ws([Ws.file] == ov.rep);
fig = wholeStackFigure(W2, ov, 4200+ov.rep);
savePNG(fig, fullfile(ov.saveDir, sprintf('ExpressionOverlap_wholeStack_rep_%02d.png', ov.rep)));

%% Section 2c : Whole structure image statistics
fig = overlapStatsFigure(T2, ov, 4300);
savePNG(fig, fullfile(ov.saveDir, 'ExpressionOverlap_wholeStack_stats.png'));

%% Section 3-0 : Glutamate time axis - camera lag correction (run before the Section 3 family)
%  Glu_Result t_ax does not mark the start of each glutamate exposure: the two cameras' motion traces
%  show that glutamate frame k covers [t_ax(k) + lag, t_ax(k) + lag + dt], lag ~ -12 ms (t_ax marks
%  the exposure end) drifting by ~ +4 ms over 10 min. Per recording:
%   - rigid shift of each camera = mean over the motion-correction patches (Volt_Result.mc,
%     Glu_Result.mc: [x patches, y patches]); the voltage shift is rotated/flipped into the glutamate
%     camera axes with the linear part of tformReg (Section 0)
%   - for each lag tau the voltage shift averaged over [t_ax + tau, t_ax + tau + dt] is correlated with
%     the glutamate shift (both minus a lagP.hp_s moving mean), per glutamate tile
%   - lag(t) = line fitted over the tiles with a correlation peak >= lagP.minR within 4 ms of their
%     median (fewer than 2: the whole-recording lag if its peak >= lagP.minR; else no correction)
%  gluLag (one entry per recording: fpath, p with lag_ms = polyval(p, t_s), tile lags) is saved to
%  GluTimeLag.mat in lagP.saveDir and used by Sections 3-3d (sta.lagCorrect): every glutamate frame
%  is placed at t_ax + lag(t) + dt/2 (exposure midpoint). Figure 7000: correlation vs lag per recording.
lagP.foi    = 14:29;
lagP.A      = tformReg.T(1:2,1:2)';   % voltage shift -> glutamate camera axes
lagP.tau_ms = -60:0.5:60;         % lags tested
lagP.hp_s   = 2;                  % motion high-pass (s)
lagP.minR   = 0.8;                % tiles used for the fit
lagP.saveDir = '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_V2+Glu/SpikeTrigGlu';
if iswindow, lagP.saveDir = macToWindowsPath(lagP.saveDir); end
if ~isfolder(lagP.saveDir), mkdir(lagP.saveDir); end
gluLag = struct('file', {}, 'fpath', {}, 'p', {}, 'info', {});
fig = figure(7000); clf; set(fig, 'Name', 'voltage-glutamate camera lag');
try, theme(fig, 'light'); catch, set(fig, 'Color', 'w'); end
tlL = tiledlayout(fig, 'flow', 'TileSpacing', 'compact', 'Padding', 'compact');
for f = lagP.foi
    fp = fpath{f};
    if ~isfile(fullfile(fp, 'Glu_Result.mat')) || ~isfile(fullfile(fp, 'Volt_Result.mat')), continue; end
    [pL, info] = gluLagFit(fp, lagP);
    gluLag(end+1) = struct('file', f, 'fpath', fp, 'p', pL, 'info', info); %#ok<SAGROW>
    fprintf('#%d glutamate lag %+.1f -> %+.1f ms (%s; whole-recording peak %+.1f ms, r %.2f)\n', f, ...
        polyval(pL, info.tileMid_s(1)), polyval(pL, info.tileMid_s(end)), info.how, info.lag_ms, info.rPeak);
    ax = nexttile(tlL); hold(ax, 'on');
    plot(ax, info.tileMid_s, info.tileLag_ms, 'o', 'Color', [0.6 0.6 0.6]);
    plot(ax, info.tileMid_s(info.good), info.tileLag_ms(info.good), 'ko', 'MarkerFaceColor', 'k');
    plot(ax, info.tileMid_s, polyval(pL, info.tileMid_s), 'r-', 'LineWidth', 1.5); hold(ax, 'off');
    ylim(ax, [-20 0]); xlabel(ax, 'time (s)'); ylabel(ax, 'lag (ms)'); box(ax, 'off');
    title(ax, sprintf('#%d (filled: used)', f), 'FontSize', 9);
end
save(fullfile(lagP.saveDir, 'GluTimeLag.mat'), 'gluLag', 'lagP');
savePNG(fig, fullfile(lagP.saveDir, 'GluTimeLag.png'));

%% Section 3 : Spike-triggered glutamate (good voltage neurons)
%  For every good neuron (Volt_Result.GoodNeuron) the glutamate of every synapse (Glu_Result dFF_glu,
%  z-scored per trace: sta.zMode) is averaged around its spikes (Volt_Result.spike). The glutamate
%  (~83 Hz, 12 ms frames, each at its exposure midpoint: t_ax + lag(t) + dt/2 with the camera lag
%  of Section 3-0, sta.lagCorrect) is linearly interpolated
%  at the voltage frame rate (~1 kHz) around every spike; no binning.
%   score : mean z of the interpolated samples inside sta.score_ms around each spike, averaged over spikes
%   p     : circular permutation test - the glutamate trace is rolled against the spikes by every
%           shift of at least sta.minShift_s (all of them, via FFT); right-tailed:
%           p = (1 + #shifts with score >= observed) / (1 + #shifts)
%   amp   : score - mean of the shifted scores (glutamate above chance before the spike)
%  Place tuning of every synapse (gluSpatialInfo, the method of si_pvalue / spatialInfo_mat from
%  Analysis_20260715 S5): glutamate events (find_spike_bh on the z-scored trace), spatial information
%  (bits/event) over running frames, and a circular-shuffle p-value.
%  Figures (sta.doPlot; rows: all spikes / isolated spikes / complex-spike onsets) show the synapses of sta.synSelect: spike-tuned (spike-triggered p < sta.alpha),
%  place-tuned (SI p < sta.siAlpha), both, or all; kymograph rows sorted by sta.synSort: x position of
%  the footprint from camera right (basal dendrite) to left, or p. One figure per neuron: kymograph,
%  mean traces (selected vs all synapses), spike-triggered voltage, and the footprints S_glu on the
%  glutamate image coloured by amp (selected synapses) and by -log10 p (all; blue circles = selected);
%  the voltage neuron's footprint is outlined in orange.
%  All values of all synapses are saved: SpikeTrigGlu.mat in sta.saveDir.
sta.foi        = [14:29];
sta.score_ms   = [-30 15];        % glutamate window around the spike for the test (ms)
sta.kymo_ms    = [-500 500];      % kymograph window (ms)
sta.base_ms    = [];     % kymograph rows are shown minus their mean here (must lie inside sta.kymo_ms;
                                  %   otherwise the samples before the test window are used, with a warning);
                                  %   [] = no baseline subtraction
sta.minShift_s = 1;               % circular shifts closer than this to zero are not used as null
sta.frameMid   = true;            % glutamate frame time = exposure midpoint (exposure start + dt/2)
sta.lagCorrect = false;           % true: glutamate exposure start = t_ax + lag(t) from Section 3-0 (gluLag)
if ~exist('gluLag', 'var')        % Section 3-0 not run in this session: its saved result
    gl = fullfile('/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_V2+Glu/SpikeTrigGlu', 'GluTimeLag.mat');
    if iswindow, gl = macToWindowsPath(gl); end
    if isfile(gl), S = load(gl); gluLag = S.gluLag; lagP = S.lagP; else, gluLag = []; end
end
sta.gluLag = gluLag;
if exist('lagP', 'var'), sta.lagP = lagP; end
sta.alpha      = 0.05;
sta.zMode      = 'robust';        % glutamate per trace: 'robust' z = (x - median) / (1.4826 MAD), 'std' = zscore
sta.silent_ms  = 30;              % isolated spikes: no other spike within +- this (ms); also Section 3c
sta.burstISI_ms = 20;             % Section 3c: burst = next spike within this (ms)
sta.pixel_um   = 0.887;           % um per glutamate-camera pixel (kymograph y axis in um from the soma)
sta.kymoClim   = [];   % kymograph colour limits (glutamate z): [] = from sta.kymoPct of each kymograph;
                                  %   [lo hi] = same for all; Section 3 figure: 3 x 2 = [all; isolated; complex]
sta.kymoPct    = [1 99];          % percentiles of the kymograph values used as colour limits when sta.kymoClim = []
sta.ampClim    = [];              % amplitude map colour limits [lo hi] (z above shuffle); [] = 20-90th pct
sta.synSelect  = 'spike';         % synapses in the figures: 'spike' (spike-triggered p < alpha), 'place', 'both', 'all'
sta.synSort    = 'x';             % kymograph rows: 'x' (camera right = basal at the top, to the left) or 'p'
sta.doPlot     = true;            % draw a figure per neuron
sta.savePlot   = true;            %   and save it as PNG in sta.saveDir
sta.siAlpha    = 0.05;
sta.siShuffle  = 1000;            % circular shuffles for the SI p-value
sta.siEventThr = [3 1];           % glutamate events: find_spike_bh [height prominence] on zscore(trace)
sta.placeBins  = 150;             % position bins over the track (as Analysis_20260715 S5)
sta.trackLength = 115;            % VR track length (VR_Result row 5 units)
sta.velThr     = 0.004;           % running: VR velocity (last row, track units per voltage frame) above this
sta.siLegacy   = false;           % true: exactly spatialInfo_mat (mean event rate over ALL frames, see gluSpatialInfo)
sta.saveDir    = '/Volumes/cohen_lab/Lab/Labmembers/Byung Hun Lee/Data/Statistics_V2+Glu/SpikeTrigGlu';
if iswindow, sta.saveDir = macToWindowsPath(sta.saveDir); end
if ~isfolder(sta.saveDir), mkdir(sta.saveDir); end

STAglu = {};
for f = sta.foi
    V = load(fullfile(fpath{f}, 'Volt_Result.mat'), 'Result'); V = V.Result;
    good = find(V.GoodNeuron(:))';
    if isempty(good), fprintf('#%d no good neurons, skipped\n', f); continue; end
    tLoad = tic;
    G = gluForSTA(fullfile(fpath{f}, 'Glu_Result.mat'), sta);
    G.dtV = median(diff(V.t_ax));                        % voltage frame interval: interpolation step
    fprintf('#%d glutamate loaded: %d synapses x %d frames (%.0f s)\n', f, G.nSyn, G.nT, toc(tLoad));
    G = addPlaceTuning(G, V, fullfile(fpath{f}, 'VR_Result.mat'), sta);
    fprintf('#%d place-tuned synapses (SI p < %g): %d / %d (%.1f expected)\n', f, sta.siAlpha, nnz(G.SIp < sta.siAlpha), ...
        G.nSyn, sta.siAlpha*G.nSyn);
    for n = good
        spk = find(V.spike(n,:) > 0);
        R3 = spikeTrigGlu(G, V.t_ax(spk), sta.score_ms, sta);
        R3.file = f; R3.roi = n; R3.mouse = raw{f,3};
        R3.SI = G.SI; R3.SIp = G.SIp; R3.nEvents = G.nEv;          % place tuning of every synapse
        STAglu{end+1} = R3; %#ok<SAGROW>
        fprintf('#%d ROI %d: %d spikes (%d used), spike-tuned %d / %d synapses (%.1f expected), of them place-tuned %d\n', ...
            f, n, numel(spk), R3.nSpk, nnz(R3.p < sta.alpha), G.nSyn, sta.alpha*G.nSyn, nnz(R3.p < sta.alpha & G.SIp < sta.siAlpha));
        if sta.doPlot
            fig = spikeTrigGluFigure(R3, G, V, n, sta, 5000 + 100*f + n, ...
                sprintf('#%d ROI %d: all spikes (%d), window [%g %g] ms', f, n, R3.nSpk, sta.score_ms));
            if sta.savePlot, savePNG(fig, fullfile(sta.saveDir, sprintf('STAglu_f%02d_roi%02d.png', f, n))); end
        end
    end
    clear G
end
STAglu = [STAglu{:}];
save(fullfile(sta.saveDir, 'SpikeTrigGlu.mat'), 'STAglu', 'sta', '-v7.3');

%% Section 3a0 : Representative recording - images and traces
%  Top: voltage reference image in the glutamate frame (Volt.ref_im_glu, cropped to the voltage camera
%  FOV) with the ROI of neuron rep.n (orange), and the glutamate mean image with the footprints of the
%  rep.nTop most spike-tuned synapses of that neuron (lowest spike-triggered p, Section 3 test with all
%  its spikes, among p < sta.alpha). Then the voltage trace of rep.n (robust z, drift removed) with its
%  detected spikes (Volt.spike, red) and complex spikes (Volt.CStrace, shaded), a kymograph of those
%  synapses' glutamate (robust z: (x - median) / 1.4826 MAD per trace, turbo; rows sorted by footprint x, camera right = basal at the top, to apical; the
%  bullet at each row has the colour of its footprint), and the VR track position (reward times marked). Traces show rep.win_s
%  ([] = whole recording); time axes linked. Uses the Section 3 parameters (sta).
rep.f        = 24;
rep.n        = 7;                 % voltage neuron shown
rep.win_s    = [];                % time window of the traces (s); [] = whole recording
rep.nTop     = 50;                % most spike-tuned synapses shown (lowest p)
rep.clim     = [-0.5 2];                % kymograph colour limits [lo hi] (robust z); [] = 0 to 99.5th pct of z
rep.markerSize = 7;               % circles marking the synapses on the glutamate image
rep.saveDir  = sta.saveDir;
if ~(exist('G', 'var') && isfield(G, 'file') && G.file == rep.f && isfield(G, 'F'))   % reuse loaded data of the same file
    V = load(fullfile(fpath{rep.f}, 'Volt_Result.mat'), 'Result'); V = V.Result;
    G = gluForSTA(fullfile(fpath{rep.f}, 'Glu_Result.mat'), sta);
    G.dtV = median(diff(V.t_ax));
    G = addPlaceTuning(G, V, fullfile(fpath{rep.f}, 'VR_Result.mat'), sta);
    G.file = rep.f;
end
Rr = spikeTrigGlu(G, V.t_ax(V.spike(rep.n,:) > 0), sta.score_ms, sta);   % spike-tuning of every synapse
tunedR = find(Rr.p < sta.alpha);
[~, o] = sortrows([Rr.p(tunedR), -Rr.amp(tunedR)]);    % lowest p first; ties (p at its floor) by amplitude
topR = tunedR(o(1:min(rep.nTop, numel(tunedR))));
[~, ox] = sort(G.xy(topR, 1), 'descend'); topR = topR(ox);   % rows: camera right (basal) at the top -> apical
VRr = load(fullfile(fpath{rep.f}, 'VR_Result.mat'), 'VR_Result'); VRr = VRr.VR_Result;
[~, recName] = fileparts(fpath{rep.f});
fig = representativeFigure(V, G, VRr, rep, topR, sprintf('%d most spike-tuned synapses (of %d, p < %g)', ...
    numel(topR), numel(tunedR), sta.alpha), 3000 + rep.f, sprintf('#%d %s, ROI %d', rep.f, recName, rep.n));
if sta.savePlot, savePNG(fig, fullfile(rep.saveDir, sprintf('Representative_f%02d_roi%02d.png', rep.f, rep.n))); end

%% Section 3a : Interactive spike-triggered glutamate (one neuron)
%  Spike type (all / simple / complex-spike onsets, Volt_Result.SpClass), minimum pre- and
%  post-ISI (ms; isolation) and the glutamate window can be changed; Update reruns the selection,
%  the kymograph and the permutation test. Pre-ISI: time since the previous spike. Post-ISI: time to
%  the next spike; for a complex-spike onset, from the end of its burst (Volt_Result.CStrace).
f = 16; n = 4;
V = load(fullfile(fpath{f}, 'Volt_Result.mat'), 'Result'); V = V.Result;
G = gluForSTA(fullfile(fpath{f}, 'Glu_Result.mat'), sta);
G.dtV = median(diff(V.t_ax));
G = addPlaceTuning(G, V, fullfile(fpath{f}, 'VR_Result.mat'), sta);
G.file = f;
interactiveSpikeTrigGlu(V, n, G, sta, sprintf('#%d ROI %d', f, n));

%% Section 3b : Event viewer - single complex / simple spikes with glutamate kymographs
%  One event at a time (left/right arrow keys or Prev/Next): top, the voltage around it (spikes red);
%  middle and bottom, kymographs of the z-scored glutamate of the "tuned" and the other synapses (at
%  the glutamate frames; rows sorted by footprint x, camera right = basal at the top, labelled in um
%  from the soma; colour limits = the set percentiles of each kymograph). Tuned = place-tuned (SI p < sta.siAlpha), spike-triggered
%  (Section 3 test, p < sta.alpha, for the chosen events) or both. "Average of all events" shows
%  the spike-triggered average (glutamate interpolated at the voltage rate) instead.
%  Isolation: only events with pre-ISI (time since the previous spike) and post-ISI (time to the next
%  spike; for a complex-spike onset, from the end of its burst) at least the set values are used, for
%  the single events, the average and the spike-triggered test.
f = 16; n = 6;
if ~(exist('G', 'var') && isfield(G, 'file') && G.file == f)   % reuse Section 3a's data when it is the same file
    V = load(fullfile(fpath{f}, 'Volt_Result.mat'), 'Result'); V = V.Result;
    G = gluForSTA(fullfile(fpath{f}, 'Glu_Result.mat'), sta);
    G.dtV = median(diff(V.t_ax));
    G = addPlaceTuning(G, V, fullfile(fpath{f}, 'VR_Result.mat'), sta);
    G.file = f;
end
interactiveEventViewer(V, n, G, sta, sprintf('#%d ROI %d', f, n));

%% Section 3c : Spike-triggered kymographs by spike class (one neuron)
%  Silent window w = sta.silent_ms:
%    i)   isolated spikes        : no other spike within +-w
%    ii)  first spike of a burst : no spike in the preceding w, next spike within sta.burstISI_ms;
%                                  spikes of complex spikes (CStrace) are excluded (they are iii)
%    iii) first spike of a complex spike : complex-spike onset (SpClass(2)) with no spike in the preceding w
%  Rows: the neuron's spike-tuned synapses (all spikes, Section 3 test, p < sta.alpha), the same rows in
%  all three panels, sorted by footprint x from camera right (basal, top) to left (apical, bottom);
%  labelled by x distance from the soma (voltage footprint centroid) in um (+ = basal side).
%  Glutamate interpolated at the voltage rate, minus the sta.base_ms baseline; one colour scale.
sta.kymoClim = [-0.3 2];
f = 24; n = 7;
if ~(exist('G', 'var') && isfield(G, 'file') && G.file == f)   % reuse the data of Section 3a/3b when it is the same file
    V = load(fullfile(fpath{f}, 'Volt_Result.mat'), 'Result'); V = V.Result;
    G = gluForSTA(fullfile(fpath{f}, 'Glu_Result.mat'), sta);
    G.dtV = median(diff(V.t_ax));
    G = addPlaceTuning(G, V, fullfile(fpath{f}, 'VR_Result.mat'), sta);
    G.file = f;
end
fig = spikeClassKymographs(V, n, G, sta, 5900 + 100*f + n, sprintf('#%d ROI %d', f, n));
if sta.savePlot, savePNG(fig, fullfile(sta.saveDir, sprintf('STAglu_classes_f%02d_roi%02d.png', f, n))); end

%% Section 3d : Spike-triggered average movies by spike class (TriggerSnippets)
%  Events of one neuron in three classes (silent window w = sta.silent_ms):
%    i)   isolated spikes  : no other spike within +-w
%    ii)  first spike of a burst : no spike in the preceding w, next spike within sta.burstISI_ms;
%                          spikes of complex spikes excluded (they are iii)
%    iii) first spike of a complex spike : complex-spike onsets (Volt.SpClass(2))
%  Their voltage and glutamate movie snippets (<fp>/TriggerSnippets, Section 6 of
%  Analysis_20261002_SD_VoltronST_plus_Glu.m; run it with opt.trig.update to add events of re-curated
%  spikes) are read one event at a time and averaged; glutamate frames (exposure midpoint) are
%  linearly interpolated onto the frame grid relative to the spike (each event has its own sub-frame
%  offset). Each pixel's average is then taken relative to its own median over time (dF = mean F -
%  median_t(mean F), no division), so static brightness differences between pixels drop out.
%  Both movies are smoothed with imgaussfilt (mov.sigma px). Voltage is shown as -dF (Voltron dims
%  with depolarisation). Viewer and .mp4 (dark): rows = classes; voltage movie, its ROI trace with a
%  moving time cursor, glutamate movie (gray). STAmovie_f<ff>_roi<nn>.mat and .mp4 in sta.saveDir.
mov.f = 16; mov.n = 4;
mov.sigma    = 1;                 % imgaussfilt of every frame (px); 0 = none
mov.pctV     = [1 99.5];          % colour limits: percentiles of each channel over all classes and frames
mov.pctG     = [1 99.5];
mov.writeMP4 = true;              % movie of the three classes (voltage frame times, every mov.step_ms)
mov.step_ms  = 1;
mov.maxEvents = Inf;              % at most this many events per class (evenly spread; Inf = all; e.g. 20 for a quick look)
SM = trigClassMovies(fpath{mov.f}, V_forMovie(fpath{mov.f}), mov.n, sta, mov);
save(fullfile(sta.saveDir, sprintf('STAmovie_f%02d_roi%02d.mat', mov.f, mov.n)), 'SM', 'mov', '-v7.3');
figM = trigMovieViewer(SM, mov, sprintf('#%d ROI %d', mov.f, mov.n));
if mov.writeMP4, trigMovieMP4(SM, mov, fullfile(sta.saveDir, sprintf('STAmovie_f%02d_roi%02d.mp4', mov.f, mov.n))); end

%% Local functions
function R = expressionOverlap(V, St, imgSize, ov)
% Somata of the voltage ROIs (orange) and of iGluSnFR (blue) in the overlapping field of one
% recording, merged into neurons and measured in both channels. Pixel indices are in the
% glutamate image (imgSize), the frame of Volt.ftprnt_glu.
H = imgSize(1); W = imgSize(2);
zc = find(St.cropSlices == St.slice);
zz = max(1, zc-ov.zHalf):min(numel(St.cropSlices), zc+ov.zHalf);
B = fovImage(max(St.cropBlue(:,:,zz), [], 3), St, H, W);
O = fovImage(max(St.cropOrange(:,:,zz), [], 3), St, H, W);
[V.ref_im_glu, V.ftprnt_glu, refReg] = registerRefToStack(V.ref_im_glu, V.ftprnt_glu, O, ov);
field = imerode(V.ref_im_glu > 0, strel('disk', ov.edge_px)) & ~isnan(B) & ~isnan(O);
Bh = highpass(B, ov.hpSigma); Oh = highpass(O, ov.hpSigma);

blue = findSomata(Bh, field, field, ov);
nV = size(V.ftprnt_glu, 3);
orange = struct('pix', cell(1,nV), 'method', '', 'ftPix', []);
for n = 1:nV
    [orange(n).pix, orange(n).method, orange(n).ftPix] = orangeSoma(Oh, V.ftprnt_glu(:,:,n), field, ov, 'orange soma', false);
end
% Somata not seen in the stack: voltage reference image in structure-orange brightness
inStack = strcmp({orange.method}, 'orange soma');
[Rn, refCal] = calibrateRef(double(V.ref_im_glu), O, vertcat(orange(inStack).pix), field, ov);
Rh = highpass(Rn, ov.hpSigma);
ovR = ov; ovR.seedZ = ov.refSeedZ;
for n = find(strcmp({orange.method}, 'none'))
    [orange(n).pix, orange(n).method] = orangeSoma(Rh, V.ftprnt_glu(:,:,n), field, ovR, 'voltage ref image', true);
end

% Same neuron: overlap / smaller area >= matchFrac, paired one-to-one from the largest overlap
M = zeros(nV, numel(blue));
for n = 1:nV
    for j = 1:numel(blue)
        a = min(numel(orange(n).pix), numel(blue(j).pix));
        if a > 0, M(n,j) = numel(intersect(orange(n).pix, blue(j).pix)) / a; end
    end
end
pairOf = zeros(1, numel(blue));
while ~isempty(M) && max(M(:)) >= ov.matchFrac
    [~, k] = max(M(:)); [n, j] = ind2sub(size(M), k);
    pairOf(j) = n; M(n,:) = 0; M(:,j) = 0;
end

% Neurons: blue somata (with their voltage ROI if matched), then unmatched voltage ROIs
nrn = struct('pix', {}, 'voltROI', {}, 'blueCell', {}, 'source', {});
for j = 1:numel(blue)
    src = 'blue'; if pairOf(j) > 0, src = 'blue+volt'; end
    nrn(end+1) = struct('pix', blue(j).pix, 'voltROI', pairOf(j), 'blueCell', j, 'source', src); %#ok<AGROW>
end
for n = setdiff(1:nV, pairOf)
    if isempty(orange(n).pix), continue; end            % footprint outside the overlapping field
    nrn(end+1) = struct('pix', orange(n).pix, 'voltROI', n, 'blueCell', 0, 'source', 'volt'); %#ok<AGROW>
end
allMask = false(H, W);
for i = 1:numel(nrn), allMask(nrn(i).pix) = true; end
for i = 1:numel(nrn)
    [r, c] = ind2sub([H W], nrn(i).pix);
    nrn(i).x = mean(c); nrn(i).y = mean(r); nrn(i).area = numel(nrn(i).pix);
    nrn(i).sRow = St.rows(1) - 1 + nrn(i).y;             % position in the whole structure stack
    nrn(i).sCol = St.cols(1) - 1 + nrn(i).x;
    nrn(i).slice = St.slice;
    nrn(i).inFOV = true;
    nrn(i).orangeMethod = '';                            % how the orange soma of a voltage ROI was found
    if nrn(i).voltROI > 0, nrn(i).orangeMethod = orange(nrn(i).voltROI).method; end
    Oimg = O; nrn(i).orangeImage = 'structure';
    if strcmp(nrn(i).orangeMethod, 'voltage ref image'), Oimg = Rn; nrn(i).orangeImage = 'voltage ref (matched)'; end
    [nrn(i).zO, nrn(i).FinO, nrn(i).FbgO] = ringContrast(Oimg, nrn(i).pix, allMask, field, ov.bgRing);
    [nrn(i).zB, nrn(i).FinB, nrn(i).FbgB] = ringContrast(B, nrn(i).pix, allMask, field, ov.bgRing);
    nrn(i).V2  = nrn(i).zO >= ov.posZ;
    nrn(i).Glu = nrn(i).zB >= ov.posZ;
end
R = struct('B', single(B), 'O', single(O), 'Ref', single(Rn), 'refCal', refCal, 'refReg', refReg, 'field', field, ...
    'blue', blue, 'orange', orange, 'neurons', nrn, 'slices', St.cropSlices(zz), 'nVoltROI', nV);
end

function W2 = wholeStackOverlap(St, R1, ov)
% Blue and orange somata detected over the whole stack slab (St.slice +-ov.zHalf), merged into
% neurons and measured in both channels. Pixel indices are in the stack image.
info = imfinfo(St.blueFile); nZ = numel(info);
zz = max(1, St.slice-ov.zHalf):min(nZ, St.slice+ov.zHalf);
B = slabMax(St.blueFile, zz); O = slabMax(St.orangeFile, zz);
sz = size(B);
valid = imerode(B > 0 & O > 0, strel('disk', ov.edge_px + 2));   % stitched canvas: zero outside the tiles
Bn = double(B); On = double(O); Bn(~valid) = NaN; On(~valid) = NaN;
ovW = ov; ovW.maxAspect = ov.wsMaxAspect;
ovW.seedZ = ov.wsSeedZ(1); blue   = findSomata(highpass(Bn, ov.hpSigma), valid, valid, ovW);
ovW.seedZ = ov.wsSeedZ(2); orange = findSomata(highpass(On, ov.hpSigma), valid, valid, ovW);

% Same neuron: overlap / smaller area >= ov.matchFrac, paired one-to-one from the largest overlap
Lb = zeros(sz); for j = 1:numel(blue), Lb(blue(j).pix) = j; end
M = zeros(numel(orange), numel(blue));
for n = 1:numel(orange)
    l = Lb(orange(n).pix); l = l(l > 0);
    for j = unique(l)'
        M(n,j) = nnz(l == j) / min(numel(orange(n).pix), numel(blue(j).pix));
    end
end
pairOf = zeros(1, numel(blue));
while ~isempty(M) && max(M(:)) >= ov.matchFrac
    [~, k] = max(M(:)); [n, j] = ind2sub(size(M), k);
    pairOf(j) = n; M(n,:) = 0; M(:,j) = 0;
end
nrn = struct('pix', {}, 'blueCell', {}, 'orangeCell', {}, 'source', {});
for j = 1:numel(blue)
    src = 'blue'; if pairOf(j) > 0, src = 'blue+orange'; end
    nrn(end+1) = struct('pix', blue(j).pix, 'blueCell', j, 'orangeCell', pairOf(j), 'source', src); %#ok<AGROW>
end
for n = setdiff(1:numel(orange), pairOf)
    nrn(end+1) = struct('pix', orange(n).pix, 'blueCell', 0, 'orangeCell', n, 'source', 'orange'); %#ok<AGROW>
end

% Voltage-recorded neurons of Section 1 (same stack coordinates)
v1 = R1.neurons([R1.neurons.voltROI] > 0);
allMask = false(sz);
for i = 1:numel(nrn), allMask(nrn(i).pix) = true; end
for i = 1:numel(nrn)
    [r, c] = ind2sub(sz, nrn(i).pix);
    nrn(i).x = mean(c); nrn(i).y = mean(r); nrn(i).area = numel(nrn(i).pix);
    nrn(i).sRow = nrn(i).y; nrn(i).sCol = nrn(i).x; nrn(i).slice = St.slice;
    nrn(i).inFOV = nrn(i).y >= St.rows(1) && nrn(i).y <= St.rows(end) && nrn(i).x >= St.cols(1) && nrn(i).x <= St.cols(end);
    nrn(i).voltROI = 0;
    if ~isempty(v1)
        [d, k] = min(hypot([v1.sRow] - nrn(i).sRow, [v1.sCol] - nrn(i).sCol));
        if d <= ov.dedupDist, nrn(i).voltROI = v1(k).voltROI; end
    end
    nrn(i).orangeMethod = 'stack detection'; nrn(i).orangeImage = 'structure';
    [nrn(i).zO, nrn(i).FinO, nrn(i).FbgO] = ringContrast(On, nrn(i).pix, allMask, valid, ov.bgRing);
    [nrn(i).zB, nrn(i).FinB, nrn(i).FbgB] = ringContrast(Bn, nrn(i).pix, allMask, valid, ov.bgRing);
    nrn(i).V2  = nrn(i).zO >= ov.posZ;
    nrn(i).Glu = nrn(i).zB >= ov.posZ;
end
W2 = struct('B', B, 'O', O, 'valid', valid, 'blue', blue, 'orangeCells', orange, 'neurons', nrn, ...
    'slices', zz, 'fovRows', St.rows([1 end]), 'fovCols', St.cols([1 end]));
end

function X = slabMax(f, zz)
% Max projection of slices zz of a tiff stack
info = imfinfo(f);
X = imread(f, zz(1), 'Info', info);
for z = zz(2:end), X = max(X, imread(f, z, 'Info', info)); end
end

function fig = wholeStackFigure(W2, ov, figNo)
% Whole stack slab: V2-ST (magenta) / iGluSnFR (green) merge, and the orange and blue channels
% with the neuron ROIs coloured by category; the recording FOV is outlined
C = categoryColors; cats = ["V2+ Glu+","V2+ only","Glu+ only","neither"];
B = double(W2.B); O = double(W2.O);
gB = mat2gray(B, prctile(B(W2.valid), [1 99.7])); gO = mat2gray(O, prctile(O(W2.valid), [1 99.7]));
fig = figure(figNo); clf; set(fig, 'Name', sprintf('#%d whole stack overlap', W2.file));
try, theme(fig, 'light'); catch, set(fig, 'Color', 'w'); end
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
pan = {cat(3, gO, gB, gO), repmat(gO, [1 1 3]), repmat(gB, [1 1 3])};
ttl = {'merge: V2-ST (magenta) / iGluSnFR (green)', 'V2-ST (orange)', 'iGluSnFR (blue)'};
fov = [W2.fovCols(1) W2.fovRows(1) diff(W2.fovCols) diff(W2.fovRows)];
for p = 1:3
    ax = nexttile(tl); image(ax, pan{p}); axis(ax, 'image', 'off'); hold(ax, 'on');
    rectangle(ax, 'Position', fov, 'EdgeColor', 'w', 'LineStyle', '--');
    if p > 1
        cn = arrayfun(@categoryOf, W2.neurons);
        for k = 1:4                                      % outlines of each category at once
            m = false(size(B)); m(vertcat(W2.neurons(cn == cats(k)).pix)) = true;
            bd = bwboundaries(m);
            for b = 1:numel(bd), plot(ax, bd{b}(:,2), bd{b}(:,1), 'Color', C(k,:), 'LineWidth', 1); end
        end
        for n = W2.neurons([W2.neurons.voltROI] > 0)
            text(ax, n.x, n.y - 14, sprintf('V%d', n.voltROI), 'Color', 'w', 'FontSize', 7, 'HorizontalAlignment', 'center');
        end
    end
    if p == 1, scaleBar(ax, size(B), ov.pixel_um, 100); end
    title(ax, ttl{p}, 'FontSize', 9);
end
h = gobjects(1, 4);
for k = 1:4, h(k) = plot(ax, NaN, NaN, '-', 'Color', C(k,:), 'LineWidth', 2); end
legend(ax, h, cats, 'Location', 'southoutside', 'NumColumns', 4, 'FontSize', 8);
title(tl, sprintf('#%d %s, slices %d-%d (dashed: recording FOV), %d neurons', W2.file, W2.structure, ...
    W2.slices([1 end]), numel(W2.neurons)), 'Interpreter', 'none', 'FontSize', 10);
drawnow;
end

function X = fovImage(I, St, H, W)
% Structure crop -> recording FOV (glutamate camera pixels); NaN where the stack does not reach
r = St.cropFOV(1) + (0:H-1); c = St.cropFOV(2) + (0:W-1);
X = nan(H, W); vr = r >= 1 & r <= size(I,1); vc = c >= 1 & c <= size(I,2);
X(vr, vc) = double(I(r(vr), c(vc)));
end

function Y = highpass(X, sig)
% Hot-pixel median filter, then the broad background removed; NaN -> 0
m = isnan(X); X(m) = median(X(~m));
X = medfilt2(X, [3 3], 'symmetric');
Y = X - imgaussfilt(X, sig); Y(m) = 0;
end

function cells = findSomata(Xh, seedMask, noiseMask, ov)
% Soma-sized bright blobs of the high-passed image Xh with their seed inside seedMask.
% Seeds: local maxima of the Hessian blobness (both curvatures negative, so dendrites are not
% seeds) whose smoothed value exceeds ov.seedZ robust SD over noiseMask; weaker seeds within
% ov.seedMinDist of a stronger one dropped. Seeds on one soma (the profile between them does not
% dip below ov.mergeFrac of the weaker peak, e.g. across the nucleus, and their peaks are within
% ov.mergeRatio) are joined into one marker.
% Each soma: its watershed basin (image smoothed by ov.wsSigma) above the local background (median
% of a ring just outside ov.maxR) + min(ov.edgeFrac * (peak - background), ov.edgeZ * ring SD),
% within ov.maxR of its seeds, opened with an ov.openR disk (strips attached dendrites), holes
% filled. When that is larger than ov.area(2) or less solid than ov.minSolidity (bright dendrites
% or a neighbour attached) or longer than ov.maxAspect x its width, the level is raised in steps up
% to ov.edgeFracMax of the peak. Somata
% cut by the edge of noiseMask (the overlapping field) are dropped (ov.edgeCut).
sz = size(Xh);
Xs = imgaussfilt(Xh, ov.somaSigma);
noise = 1.4826 * mad(Xs(noiseMask), 1);
bl = hessianBlobness(Xs);
cand = find(imregionalmax(bl) & bl > 0 & Xs > ov.seedZ*noise & seedMask);
[~, o] = sort(Xs(cand), 'descend'); cand = cand(o);
[r, c] = ind2sub(sz, cand);
keep = true(size(cand));
for i = 2:numel(cand)
    d = hypot(r(1:i-1) - r(i), c(1:i-1) - c(i));
    keep(i) = ~any(d(keep(1:i-1)) < ov.seedMinDist);
end
cand = cand(keep); r = r(keep); c = c(keep);
cells = struct('pix', {}, 'seed', {}, 'peak', {});
if isempty(cand), return; end
X1 = imgaussfilt(Xh, 1);
pad = ov.maxR + 10;                                      % every step works in a window this far around a soma
[dc, dr] = meshgrid(-2:2); d2 = hypot(dr, dc) <= 2;

% join seeds of one soma (strongest first, so each group is headed by its brightest seed)
grp = 1:numel(cand);
pk = zeros(numel(cand), 1);
for i = 1:numel(cand)                                    % peak: max within 2 px of the seed
    rr = r(i) + dr(d2); cc = c(i) + dc(d2); ok = rr >= 1 & rr <= sz(1) & cc >= 1 & cc <= sz(2);
    pk(i) = max(X1(sub2ind(sz, rr(ok), cc(ok))));
end
for i = 2:numel(cand)
    for j = 1:i-1
        if grp(j) ~= j || hypot(r(i) - r(j), c(i) - c(j)) > 2*ov.maxR, continue; end
        t = linspace(0, 1, ceil(hypot(r(i) - r(j), c(i) - c(j))) + 1);
        prof = X1(sub2ind(sz, round(r(j) + t*(r(i) - r(j))), round(c(j) + t*(c(i) - c(j)))));
        bg = median(ringValues(X1, noiseMask, r([i j]), c([i j]), ov));
        if min(prof) - bg >= ov.mergeFrac * (min(pk([i j])) - bg) && ...
                min(pk([i j])) - bg >= ov.mergeRatio * (max(pk([i j])) - bg)
            grp(i) = j;
            break
        end
    end
end
mkG = num2cell(cand(:))';                                % marker pixels of each group (seeds + joining lines)
for i = find(grp ~= 1:numel(grp))
    j = grp(i); t = linspace(0, 1, ceil(hypot(r(i) - r(j), c(i) - c(j))) + 1);
    mkG{j} = [mkG{j}; sub2ind(sz, round(r(j) + t'*(r(i) - r(j))), round(c(j) + t'*(c(i) - c(j))))];
end
mk = false(sz); mk(vertcat(mkG{:})) = true;
L = watershed(imimposemin(-imgaussfilt(Xh, ov.wsSigma), mk));   % smoothed: no split at the nucleus
for g = find(grp == 1:numel(grp))
    [gr, gc] = ind2sub(sz, mkG{g});
    wr = max(1, min(gr)-pad):min(sz(1), max(gr)+pad); wc = max(1, min(gc)-pad):min(sz(2), max(gc)+pad);
    ms = false(numel(wr), numel(wc)); ms(sub2ind(size(ms), gr-wr(1)+1, gc-wc(1)+1)) = true;
    Lw = L(wr, wc); Xw = X1(wr, wc); Nw = noiseMask(wr, wc);
    lab = unique(Lw(ms)); lab = lab(lab > 0);
    if isempty(lab), continue; end
    v = ringValues(X1, noiseMask, gr, gc, ov); bg = median(v);
    if pk(g) <= bg, continue; end
    lev0 = bg + min(ov.edgeFrac * (pk(g) - bg), ov.edgeZ * 1.4826 * mad(v, 1));
    levs = linspace(lev0, max(lev0, bg + ov.edgeFracMax * (pk(g) - bg)), 6);
    inR = ismember(Lw, lab) & bwdist(ms) <= ov.maxR;
    r0 = r(g) - wr(1) + 1; c0 = c(g) - wc(1) + 1;
    for lev = levs                                       % lowest edge level giving a compact soma
        m = bwselect((inR & Xw >= lev) | ms, c0, r0, 8);
        m = imfill(imopen(m, strel('disk', ov.openR)) | ms, 'holes');   % strip dendrites, keep the seeds
        m = bwselect(m, c0, r0, 8);
        p = regionprops(m, 'Area', 'Solidity', 'MajorAxisLength', 'MinorAxisLength');
        if p.Area < ov.area(1), break; end
        if p.Area <= ov.area(2) && p.Solidity >= ov.minSolidity && p.MajorAxisLength <= ov.maxAspect * p.MinorAxisLength
            if nnz(imdilate(m, strel('disk', 2)) & ~Nw) > ov.edgeCut * p.Area, break; end   % cut by the field edge
            [mr, mc] = find(m);
            cells(end+1) = struct('pix', sub2ind(sz, mr+wr(1)-1, mc+wc(1)-1), 'seed', cand(g), 'peak', Xs(cand(g))); %#ok<AGROW>
            break
        end
    end
end
end

function v = ringValues(X1, noiseMask, rs, cs, ov)
% X1 inside noiseMask in the ring ov.maxR < distance <= ov.maxR + 8 around the pixels (rs, cs)
sz = size(X1); pad = ov.maxR + 10;
wr = max(1, min(rs)-pad):min(sz(1), max(rs)+pad); wc = max(1, min(cs)-pad):min(sz(2), max(cs)+pad);
m = false(numel(wr), numel(wc)); m(sub2ind(size(m), rs(:)-wr(1)+1, cs(:)-wc(1)+1)) = true;
d = bwdist(m); x = X1(wr, wc);
v = x(d > ov.maxR & d <= ov.maxR + 8 & noiseMask(wr, wc));
end

function b = hessianBlobness(Xs)
% Minus the larger Hessian eigenvalue, kept where both are negative (bright blob)
[gx, gy] = gradient(Xs);
[gxx, gxy] = gradient(gx); [~, gyy] = gradient(gy);
tr = gxx + gyy; dt = gxx.*gyy - gxy.^2;
b = max(-(tr/2 + sqrt(max(tr.^2/4 - dt, 0))), 0);
end

function [pix, method, ftPix] = orangeSoma(Oh, w, field, ov, label, allowDisk)
% Orange soma of one voltage ROI in the high-passed image Oh: the soma-sized blob inside the
% footprint region with the highest (peak x mean normalised footprint weight), method = label.
% None: a disk at the footprint peak (allowDisk) or empty with method 'none'.
ws = imgaussfilt(double(w), 2); ws(~field) = 0;
pix = []; ftPix = []; method = 'outside field';
if max(ws(:)) <= 0, return; end
[~, ip] = max(ws(:)); [r0, c0] = ind2sub(size(ws), ip);
mf = bwselect(ws > ov.ftThr * max(ws(:)), c0, r0, 8);
ftPix = find(mf);
cand = findSomata(Oh, mf, field, ov);
if ~isempty(cand)
    wn = ws / max(ws(:));
    [~, k] = max(arrayfun(@(x) x.peak * mean(wn(x.pix)), cand));
    pix = cand(k).pix; method = label;
elseif ~allowDisk
    method = 'none';
else
    [cc, rr] = meshgrid(1:size(ws,2), 1:size(ws,1));
    pix = find(hypot(rr - r0, cc - c0) <= ov.somaR & field); method = 'footprint peak disk';
end
end

function [Rv, Fp, reg] = registerRefToStack(Rv, Fp, O, ov)
% Translation of the voltage reference image onto the stack's orange channel (normalised cross-
% correlation of band-passed images, +-ov.refMaxShift px); the footprints Fp (same camera, same
% transform) are shifted with it. Not applied when the peak is below ov.refMinCorr or on the edge
% of the search window.
Rv = double(Rv);
in = imerode(Rv > 0, strel('disk', 5));                  % voltage camera FOV without its edge
[r, c] = find(in);
rr = min(r):max(r); cc = min(c):max(c);
X = Rv; X(~(Rv > 0)) = NaN;
T = imgaussfilt(highpass(X, ov.hpSigma), 1.5); T(~in) = 0; T = T(rr, cc);
A = imgaussfilt(highpass(O, ov.hpSigma), 1.5);
L = ov.refMaxShift;
ar = rr(1)-L:rr(end)+L; ac = cc(1)-L:cc(end)+L;
Ap = zeros(numel(ar), numel(ac));
vr = ar >= 1 & ar <= size(A,1); vc = ac >= 1 & ac <= size(A,2);
Ap(vr, vc) = A(ar(vr), ac(vc));
C = normxcorr2(T, Ap);
C = C(size(T,1):end-size(T,1)+1, size(T,2):end-size(T,2)+1);   % template fully inside: shifts -L..L
[pk, ix] = max(C(:)); [i, j] = ind2sub(size(C), ix);
dy = i - 1 - L; dx = j - 1 - L;
C0 = C(L+1, L+1);                                        % correlation without a shift
reg = struct('dx', dx, 'dy', dy, 'corr', pk, 'corr0', C0, 'applied', pk >= ov.refMinCorr && abs(dx) < L && abs(dy) < L);
if ~reg.applied, return; end
Rv = imtranslate(Rv, [dx dy], 'FillValues', 0);
Fp = imtranslate(double(Fp), [dx dy], 'FillValues', 0);
end

function [Rn, cal] = calibrateRef(Rv, O, brightPix, field, ov)
% Voltage reference image (Volt.ref_im_glu) in the brightness of the structure's orange channel:
% gain and offset matching the background (median over the field) and the bright cells (median
% over the orange somata found in the stack; 99.5th percentiles when they cover < ov.calMinPix px)
ok = field & Rv > 0;
bgR = median(Rv(ok)); bgO = median(O(ok));
src = 'orange somata in the stack';
if numel(brightPix) >= ov.calMinPix
    brR = median(Rv(brightPix)); brO = median(O(brightPix));
end
if numel(brightPix) < ov.calMinPix || brR <= bgR || brO <= bgO
    brR = prctile(Rv(ok), 99.5); brO = prctile(O(ok), 99.5); src = '99.5th percentiles';
end
gain = (brO - bgO) / (brR - bgR);
Rn = bgO + gain * (Rv - bgR);
Rn(~(Rv > 0)) = NaN;                                     % outside the voltage camera FOV
cal = struct('gain', gain, 'offset', bgO - gain*bgR, 'source', src);
end

function [z, Fin, Fbg] = ringContrast(I, pix, allMask, field, ring)
% Mean in the ROI against a ring ring(1)-ring(2) px outside it (other ROIs and outside the field
% excluded): z = (mean - ring median) / ring robust SD
sz = size(I); [r, c] = ind2sub(sz, pix); pad = ring(2) + 2;    % window around the ROI
wr = max(1, min(r)-pad):min(sz(1), max(r)+pad); wc = max(1, min(c)-pad):min(sz(2), max(c)+pad);
m = false(numel(wr), numel(wc)); m(sub2ind(size(m), r-wr(1)+1, c-wc(1)+1)) = true;
rg = imdilate(m, strel('disk', ring(2))) & ~imdilate(m, strel('disk', ring(1))) & ~allMask(wr,wc) & field(wr,wc);
Iw = I(wr, wc); v = Iw(rg);
Fin = mean(I(pix), 'omitnan'); Fbg = median(v);
z = (Fin - Fbg) / (1.4826 * mad(v, 1) + eps);
end

function T = overlapTable(Ov, ov)
% One row per neuron and recording; dup = same neuron as in an earlier recording of the same
% structure stack (within ov.dedupDist px and ov.dedupDz slices)
rows = {};
for k = 1:numel(Ov)
    R = Ov(k);
    for i = 1:numel(R.neurons)
        n = R.neurons(i);
        rows(end+1,:) = {R.file, R.mouse, string(R.structure), i, string(n.source), n.voltROI, n.blueCell, ...
            n.x, n.y, n.sRow, n.sCol, n.slice, n.area, n.zO, n.zB, n.V2, n.Glu, string(n.orangeMethod), ...
            string(n.orangeImage), n.inFOV}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'file','mouse','structure','neuron','source','voltROI','blueCell', ...
    'x','y','sRow','sCol','slice','area','zO','zB','V2','Glu','orangeMethod','orangeImage','inFOV'});
T.recorded = T.voltROI > 0;
T.dup = false(height(T),1);
T.recordedAny = T.recorded;                              % recorded in any copy (kept on the first copy)
for i = 2:height(T)
    j = find(T.structure(1:i-1) == T.structure(i) & T.file(1:i-1) ~= T.file(i) & ~T.dup(1:i-1) & ...
        hypot(T.sRow(1:i-1) - T.sRow(i), T.sCol(1:i-1) - T.sCol(i)) <= ov.dedupDist & ...
        abs(T.slice(1:i-1) - T.slice(i)) <= ov.dedupDz, 1);
    if ~isempty(j), T.dup(i) = true; T.recordedAny(j) = T.recordedAny(j) | T.recorded(i); end
end
lab = repmat("neither", height(T), 1);
lab(T.V2 & T.Glu) = "V2+ Glu+"; lab(T.V2 & ~T.Glu) = "V2+ only"; lab(~T.V2 & T.Glu) = "Glu+ only";
T.category = categorical(lab, ["V2+ Glu+","V2+ only","Glu+ only","neither"]);
end

function printOverlapStats(U, T, ttl)
if nargin < 3, ttl = 'overlapping field, structure stack'; end
fprintf('\n=== V2-ST / iGluSnFR expression overlap (%s) ===\n', ttl);
fprintf('%d neuron entries in %d recordings, %d unique neurons in %d structure stacks\n', ...
    height(T), numel(unique(T.file)), height(U), numel(unique(U.structure)));
c = countcats(U.category); nm = categories(U.category);
for k = 1:numel(nm), fprintf('  %-10s %4d (%4.1f%%)\n', nm{k}, c(k), 100*c(k)/height(U)); end
fprintf('  Glu+ among V2+ neurons            : %d / %d (%.1f%%)\n', nnz(U.V2 & U.Glu), nnz(U.V2), 100*nnz(U.V2 & U.Glu)/max(1,nnz(U.V2)));
fprintf('  V2+ among Glu+ neurons            : %d / %d (%.1f%%)\n', nnz(U.V2 & U.Glu), nnz(U.Glu), 100*nnz(U.V2 & U.Glu)/max(1,nnz(U.Glu)));
fprintf('  Glu+ among voltage-recorded       : %d / %d (%.1f%%)\n', nnz(U.recordedAny & U.Glu), nnz(U.recordedAny), 100*nnz(U.recordedAny & U.Glu)/max(1,nnz(U.recordedAny)));
fprintf('  V2+ among voltage-recorded (QC)   : %d / %d (%.1f%%)\n', nnz(U.recordedAny & U.V2), nnz(U.recordedAny), 100*nnz(U.recordedAny & U.V2)/max(1,nnz(U.recordedAny)));
for mth = intersect(["orange soma", "voltage ref image", "footprint peak disk"], unique(U.orangeMethod)', 'stable')
    d = U.orangeMethod == mth;
    fprintf('  voltage ROI soma from %-20s: %3d, V2+ %3d, Glu+ %3d\n', mth, nnz(d), nnz(d & U.V2), nnz(d & U.Glu));
end
end

function C = categoryColors
% V2+ Glu+, V2+ only, Glu+ only, neither
C = [hex2rgb('#2a78d6'); hex2rgb('#eb6834'); hex2rgb('#1baf7a'); 0.6 0.6 0.6];
end

function c = hex2rgb(h)
c = sscanf(h(2:end), '%2x%2x%2x')' / 255;
end

function fig = overlapRepFigure(R, ov, figNo)
% Recording FOV merge (overlapping field outlined), then the overlapping field: voltage reference
% image (in structure-orange brightness, same display range as the stack's orange), orange and
% blue with the neuron ROIs coloured by category (dashed: soma found in the reference image), and
% the merge with voltage ROI numbers
C = categoryColors; cats = ["V2+ Glu+","V2+ only","Glu+ only","neither"];
B = double(R.B); O = double(R.O); Rf = double(R.Ref);
[rr, cc] = find(R.field);
yr = max(1, min(rr)-10):min(size(B,1), max(rr)+10); xr = max(1, min(cc)-10):min(size(B,2), max(cc)+10);
gB = @(X) mat2gray(X, prctile(B(R.field), [1 99.7]));
gO = @(X) mat2gray(X, prctile(O(R.field), [1 99.7]));
mergeRGB = @(o, b) cat(3, o, b, o);                      % orange -> magenta, blue -> green
fig = figure(figNo); clf; set(fig, 'Name', sprintf('#%d expression overlap', R.file));
try, theme(fig, 'light'); catch, set(fig, 'Color', 'w'); end
tl = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(tl, 1, [1 4]);
Bn = gB(B); On = gO(O); Bn(isnan(B)) = 0; On(isnan(O)) = 0;
image(ax, mergeRGB(On, Bn)); axis(ax, 'image', 'off'); hold(ax, 'on');
bd = bwboundaries(R.field);
for b = 1:numel(bd), plot(ax, bd{b}(:,2), bd{b}(:,1), 'w--', 'LineWidth', 1); end
scaleBar(ax, size(B), ov.pixel_um, 50);
title(ax, sprintf('#%d %s, slices %d-%d: V2-ST (magenta) / iGluSnFR (green); dashed = overlapping field', ...
    R.file, R.structure, R.slices([1 end])), 'Interpreter', 'none', 'FontSize', 9);
pan = {{'voltage ref image', sprintf('shift [%+d %+d] px, brightness matched', R.refReg.dx*R.refReg.applied, R.refReg.dy*R.refReg.applied)}, ...
    'V2-ST (orange, stack)', 'iGluSnFR (blue, stack)', {'merge (stack)', 'dashed: orange from ref image'}};
for p = 1:4
    ax = nexttile(tl, 4+p);
    switch p
        case 1, X = gO(Rf(yr,xr)); X(isnan(X)) = 0; image(ax, repmat(X, [1 1 3]));
        case 2, image(ax, repmat(gO(O(yr,xr)), [1 1 3]));
        case 3, image(ax, repmat(gB(B(yr,xr)), [1 1 3]));
        case 4, image(ax, mergeRGB(gO(O(yr,xr)), gB(B(yr,xr))));
    end
    axis(ax, 'image', 'off'); hold(ax, 'on');
    for i = 1:numel(R.neurons)
        n = R.neurons(i);
        ci = find(cats == categoryOf(n));
        m = false(size(B)); m(n.pix) = true; m = m(yr, xr);
        ls = '-'; if strcmp(n.orangeImage, 'voltage ref (matched)'), ls = '--'; end
        if p < 4
            bd = bwboundaries(m);
            for b = 1:numel(bd), plot(ax, bd{b}(:,2), bd{b}(:,1), ls, 'Color', C(ci,:), 'LineWidth', 1.5); end
        end
        if n.voltROI > 0
            text(ax, n.x - xr(1) + 1, n.y - yr(1) + 1 - 12, sprintf('V%d', n.voltROI), 'Color', 'w', ...
                'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end
    end
    if p == 4, scaleBar(ax, [numel(yr) numel(xr)], ov.pixel_um, 20); end
    title(ax, pan{p}, 'FontSize', 9);
end
h = gobjects(1, 4);                                      % ROI category legend under the merge
for k = 1:4, h(k) = plot(ax, NaN, NaN, '-', 'Color', C(k,:), 'LineWidth', 2); end
legend(ax, h, cats, 'Location', 'southoutside', 'NumColumns', 2, 'FontSize', 8);
drawnow;
end

function s = categoryOf(n)
if n.V2 && n.Glu, s = "V2+ Glu+"; elseif n.V2, s = "V2+ only"; elseif n.Glu, s = "Glu+ only"; else, s = "neither"; end
end

function scaleBar(ax, sz, um_px, um)
L = um / um_px; x0 = sz(2) - L - 10; y0 = sz(1) - 10;
plot(ax, [x0 x0+L], [y0 y0], 'w-', 'LineWidth', 3);
text(ax, x0 + L/2, y0 - 8, sprintf('%d \\mum', um), 'Color', 'w', 'FontSize', 8, 'HorizontalAlignment', 'center');
end

function fig = overlapStatsFigure(T, ov, figNo)
% Unique neurons: orange vs blue contrast, and the overlap fractions pooled with one dot per mouse
C = categoryColors; U = T(~T.dup, :);
fig = figure(figNo); clf; set(fig, 'Name', 'V2-ST / iGluSnFR overlap');
try, theme(fig, 'light'); catch, set(fig, 'Color', 'w'); end
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% 1: z orange vs z blue, all unique neurons alike
ax = nexttile(tl); hold(ax, 'on');
scatter(ax, U.zO, U.zB, 30, C(1,:), 'filled', 'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'w');
xline(ax, ov.posZ, ':', 'Color', [0.4 0.4 0.4]); yline(ax, ov.posZ, ':', 'Color', [0.4 0.4 0.4]);
xlabel(ax, 'V2-ST contrast z (orange)'); ylabel(ax, 'iGluSnFR contrast z (blue)');
title(ax, sprintf('%d unique neurons (dotted: z = %g)', height(U), ov.posZ), 'FontSize', 9);
box(ax, 'off');
% 2: pooled fractions, one dot per mouse
ax = nexttile(tl); hold(ax, 'on');
grp = {U.V2, U.Glu, U.recordedAny};
num = {U.Glu, U.V2, U.Glu};
lbl = ["Glu+ | V2+", "V2+ | Glu+", "Glu+ | recorded"];
mice = unique(U.mouse);
for g = 1:3
    p = 100 * nnz(grp{g} & num{g}) / max(1, nnz(grp{g}));
    bar(ax, g, p, 0.6, 'FaceColor', [0.85 0.85 0.85], 'EdgeColor', 'none');
    text(ax, g, 4, sprintf('%d/%d', nnz(grp{g} & num{g}), nnz(grp{g})), 'HorizontalAlignment', 'center', 'FontSize', 8);
    for m = 1:numel(mice)
        s = U.mouse == mice(m) & grp{g};
        if nnz(s) == 0, continue; end
        plot(ax, g + 0.15*(m - (numel(mice)+1)/2)/numel(mice), 100*mean(num{g}(s)), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    end
end
set(ax, 'XTick', 1:3, 'XTickLabel', lbl); ylim(ax, [0 110]); ylabel(ax, '% of neurons');
title(ax, {'overlap (bars: pooled unique neurons)', sprintf('dots: mice (%s)', strjoin(string(mice), ', '))}, 'FontSize', 9);
box(ax, 'off');
drawnow;
end

function savePNG(fig, fn)
try
    exportgraphics(fig, fn, 'Resolution', 150);
catch
    try, print(fig, fn, '-dpng', '-r150'); catch ME, warning('PNG export failed for %s : %s', fn, ME.message); end
end
end

function G = gluForSTA(gf, sta)
% Glutamate traces z-scored per synapse (sta.zMode; single), frame times (exposure midpoint when
% sta.frameMid), footprints as sparse [pixels x synapses], mean image, and the FFT of every trace
% along time (for the permutation test)
S = loadRetry(gf, 'dFF_glu', 't_ax', 'S_glu', 'AvgGluImg', 'frameRate');
X = single(S.dFF_glu); S = rmfield(S, 'dFF_glu');
X(~isfinite(X)) = 0;
if strcmp(sta.zMode, 'robust')
    mu = median(X, 2); sd = 1.4826 * mad(X, 1, 2);
else
    mu = mean(X, 2); sd = std(X, 0, 2);
end
G.X = (X - mu) ./ max(sd, eps('single')); clear X
[G.nSyn, G.nT] = size(G.X);
G.t = S.t_ax(1:G.nT); G.t = G.t(:)';
G.dt = median(diff(G.t));
G.lagFit = [0 0]; G.lagInfo = struct();
if isfield(sta, 'lagCorrect') && sta.lagCorrect           % exposure start = t_ax + lag(t)
    G.lagFit = gluLagOf(fileparts(gf), sta);
    G.t = G.t + polyval(G.lagFit, G.t) / 1000;
end
if sta.frameMid, G.t = G.t + G.dt/2; end
G.img = S.AvgGluImg;
Sg = S.S_glu;
if ndims(Sg) == 3 || size(Sg,1) ~= numel(G.img), Sg = reshape(Sg, numel(G.img), []); end
G.S = sparse(double(Sg));
[pix, k, w] = find(G.S); [r, c] = ind2sub(size(G.img), pix);
ws = accumarray(k, w, [size(G.S,2) 1]);
G.xy = [accumarray(k, w.*c, [size(G.S,2) 1]) accumarray(k, w.*r, [size(G.S,2) 1])] ./ max(ws, eps);   % footprint centroids [x y]
if ~isfield(sta, 'fft') || sta.fft
    G.F = fft(G.X, [], 2);                             % nSyn x nT, complex single (permutation test)
end
end

function R = spikeTrigGlu(G, tsp, score_ms, sta)
% Spike-triggered glutamate of spike times tsp (s). The glutamate traces are linearly interpolated
% at the voltage frame rate (G.dtV) around every spike - interpolation is linear in the trace, so
% this is done with sparse weights (two neighbouring glutamate frames per sample) instead of
% interpolating the whole recording:
%   STA   : kymograph, every voltage frame in sta.kymo_ms (R.t_ms)
%   score : mean of the interpolated samples in score_ms, averaged over spikes
%   p     : circular permutation over every glutamate-frame shift >= sta.minShift_s (score at shift
%           s = circular cross-correlation of the trace with the weights, all shifts at once by FFT)
tsp = tsp(:);
lo = min(score_ms(1), sta.kymo_ms(1)) / 1000; hi = max(score_ms(2), sta.kymo_ms(2)) / 1000;
used = tsp + lo >= G.t(1) & tsp + hi <= G.t(end);       % whole window inside the glutamate recording
tsp = tsp(used);
R = struct('nSpk', numel(tsp), 'score_ms', score_ms, 'used', used');
off = round(sta.kymo_ms(1)/1000/G.dtV):round(sta.kymo_ms(2)/1000/G.dtV);   % voltage frames
R.t_ms = off * G.dtV * 1000;
if isempty(tsp)
    R.STA = nan(G.nSyn, numel(R.t_ms), 'single'); [R.score, R.amp, R.p, R.nullMean, R.nullSD] = deal(nan(G.nSyn,1)); R.nNull = 0; return
end
W = interpWeights(G, tsp, off * G.dtV);                  % nT x nOff, each column averages the spikes
R.STA = zeros(G.nSyn, numel(off), 'single');
for r = 1:256:G.nSyn                                     % chunks: sparse products need double
    rr = r:min(G.nSyn, r+255);
    R.STA(rr,:) = single(double(G.X(rr,:)) * W);
end

% score (interpolated samples in score_ms) and its permutation null
offS = ceil(score_ms(1)/1000/G.dtV):floor(score_ms(2)/1000/G.dtV);
w = full(sum(interpWeights(G, tsp, offS * G.dtV), 2))' / numel(offS);
Fw = conj(fft(single(w)));
shifts = (0:G.nT-1) * G.dt;
useNull = shifts >= sta.minShift_s & shifts <= (G.nT*G.dt - sta.minShift_s);
[R.score, R.nullMean, R.nullSD, R.p] = deal(nan(G.nSyn, 1));
for r = 1:256:G.nSyn
    rr = r:min(G.nSyn, r+255);
    c = real(ifft(G.F(rr,:) .* Fw, [], 2));              % c(:,s+1) = score with the trace rolled by s
    nul = c(:, useNull);
    R.score(rr) = c(:,1);
    R.nullMean(rr) = mean(nul, 2); R.nullSD(rr) = std(nul, 0, 2);
    R.p(rr) = (1 + sum(nul >= c(:,1), 2)) / (1 + nnz(useNull));
end
R.amp = R.score - R.nullMean;
R.nNull = nnz(useNull);
end

function W = interpWeights(G, tsp, offs)
% Sparse [glutamate frames x offsets] weights of linear interpolation at times tsp + offs (s),
% averaged over the spikes: X * W(:,j) = mean over spikes of the trace interpolated at tsp + offs(j)
pos = (tsp(:) + offs(:)' - G.t(1)) / G.dt + 1;           % fractional glutamate frame (frames are evenly spaced)
k = min(max(floor(pos), 1), G.nT - 1);
a = min(max(pos - k, 0), 1);
col = repmat(1:numel(offs), numel(tsp), 1);
W = sparse([k(:); k(:) + 1], [col(:); col(:)], [1 - a(:); a(:)], G.nT, numel(offs)) / numel(tsp);
end

function fig = spikeTrigGluFigure(R, G, V, n, sta, figNo, ttl)
% Section 3 figure, three rows sharing the synapse pool and row order (sta.synSelect / sta.synSort,
% chosen with all spikes) and one kymograph colour scale:
%   top    : all spikes - kymograph, amplitude map, -log10 p map
%   middle : isolated spikes (no other spike within +-sta.silent_ms) - kymograph, mean glutamate
%            (pool vs all synapses), spike-triggered voltage
%   bottom : complex-spike onsets (Volt.SpClass(2)) - the same
% Each kymograph has its own colour limits (sta.kymoClim: [] = its sta.kymoPct percentiles, 1 x 2 shared, 3 x 2 per row).
ftV = V.ftprnt_glu(:,:,n);
fig = figure(figNo); clf; set(fig, 'Name', ttl);
try, theme(fig, 'light'); catch, set(fig, 'Color', 'w'); end
tl = tiledlayout(fig, 6, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
axK = [nexttile(tl, 1, [2 2]) nexttile(tl, 9, [2 2]) nexttile(tl, 17, [2 2])];   % kymographs (left half)
axM = [nexttile(tl, 3, [1 2]) nexttile(tl, 7, [1 2])];                          % amplitude, p maps (top right)
axG = [nexttile(tl, 11, [2 1]) nexttile(tl, 19, [2 1])]; axV = [nexttile(tl, 12, [2 1]) nexttile(tl, 20, [2 1])];
C = classifySpikes(V, n, sta.silent_ms, sta.burstISI_ms);
cls = {C.isolated, find(V.SpClass(2,:,n) > 0)};
nm = {sprintf('isolated (\\pm%g ms)', sta.silent_ms), 'complex-spike onsets'};
Rc = cell(1, 2);
for c = 1:2, Rc{c} = spikeTrigGlu(G, V.t_ax(cls{c}), sta.score_ms, sta); end
sel = []; if R.nSpk > 0, sel = synSelection(R, G, sta); end
% colour limits per kymograph (all / isolated / complex): sta.kymoClim [] = each automatic,
% [lo hi] = the same for all three, 3 x 2 = one row per kymograph
kc = sta.kymoClim; if size(kc, 1) == 1, kc = repmat(kc, 3, 1); end
s2 = sta; s2.kymoClim = []; if ~isempty(kc), s2.kymoClim = kc(1,:); end
spikeTrigGluPanels(R, G, ftV, s2, [axK(1) axM]);
title(axK(1), sprintf('all spikes (%d): %d %s synapses (%.1f expected)', R.nSpk, nnz(sel), selName(sta), sta.alpha*G.nSyn), 'FontSize', 9);
ylabel(axK(1), 'synapse (basal \rightarrow apical, \mum from soma)');
trV = voltZ(V, n); dtV = median(diff(V.t_ax)) * 1000;
for c = 1:2
    ax = axK(c+1);
    if Rc{c}.nSpk == 0 || isempty(sel) || ~any(sel)
        cla(ax); text(ax, 0.5, 0.5, sprintf('no %s', nm{c}), 'Units', 'normalized', 'HorizontalAlignment', 'center');
        cla(axG(c)); cla(axV(c)); continue
    end
    [ord, tk, tkLbl, sortLbl] = synOrder(R, G, ftV, sta, sel);
    base = staBaseline(Rc{c}, sta);
    Kc = Rc{c}.STA(ord,:) - base(ord);
    s3 = sta; s3.kymoClim = []; if ~isempty(kc), s3.kymoClim = kc(c+1,:); end
    imagesc(ax, Rc{c}.t_ms, 1:numel(ord), Kc, kymoClim(Kc, s3)); colormap(ax, kymoMap); hold(ax, 'on');
    xline(ax, 0, 'k-'); xline(ax, sta.score_ms, 'k--'); hold(ax, 'off');
    xlim(ax, Rc{c}.t_ms([1 end]));
    if strcmp(sta.synSort, 'x'), set(ax, 'YTick', tk, 'YTickLabel', tkLbl); end
    cb = colorbar(ax); cb.Label.String = gluLabel(sta);
    xlabel(ax, 'time from spike (ms)');
    if strcmp(sta.synSort, 'x'), ylabel(ax, 'synapse (basal \rightarrow apical, \mum from soma)'); else, ylabel(ax, sprintf('synapse (%s)', sortLbl)); end
    title(ax, sprintf('%s (%d spikes), same synapses', nm{c}, Rc{c}.nSpk), 'FontSize', 9);
    meanTracePanel(axG(c), Rc{c}, Rc{c}.STA - base, sta, sel, selName(sta));
    title(axG(c), sprintf('glutamate, %s', nm{c}), 'FontSize', 9);
    voltTracePanel(axV(c), trV, cls{c}(Rc{c}.used), dtV, sta, sta.score_ms);
    title(axV(c), sprintf('voltage, %s', nm{c}), 'FontSize', 9);
end
linkaxes([axK axG axV], 'x');
title(tl, ttl, 'Interpreter', 'none', 'FontSize', 10);
drawnow;
end

function tr = voltZ(V, n)
% Voltage trace of ROI n, slow drift removed (1000-frame moving median), in robust SD units
tr = double(V.traces(n,:)); tr = tr - movmedian(tr, 1000);
tr = tr / (1.4826 * mad(tr, 1));
end

function voltTracePanel(ax, trV, spk, dtV, sta, win)
% Spike-triggered voltage: up to 50 single spikes (grey) and the mean (black)
o = round(sta.kymo_ms(1)/dtV):round(sta.kymo_ms(2)/dtV);
spk = spk(spk + o(1) >= 1 & spk + o(end) <= numel(trV));
cla(ax); hold(ax, 'on');
if ~isempty(spk)
    M = trV(spk(:) + o);
    plot(ax, o*dtV, M(1:min(end,50),:)', 'Color', [0.8 0.8 0.8]);
    plot(ax, o*dtV, mean(M, 1), 'k', 'LineWidth', 1.5);
end
xline(ax, 0, 'k-'); xline(ax, win, 'k--'); hold(ax, 'off'); xlim(ax, sta.kymo_ms); box(ax, 'off');
xlabel(ax, 'time from spike (ms)'); ylabel(ax, 'voltage (robust z)');
title(ax, sprintf('spike-triggered voltage, %d spikes (grey: up to 50, black: mean)', numel(spk)), 'FontSize', 9);
end

function spikeTrigGluPanels(R, G, ftV, sta, ax)
% ax(1) kymograph, ax(2) amplitude map, ax(3) -log10 p map, ax(4) mean traces (axes or uiaxes).
% Kymograph and amplitude map: the synapses of sta.synSelect, rows sorted by sta.synSort.
sel = [];
if R.nSpk > 0, sel = synSelection(R, G, sta); end
if R.nSpk == 0 || ~any(sel)                              % nothing to show: clear and say so
    msg = 'no spikes selected'; if R.nSpk > 0, msg = sprintf('no %s synapses', selName(sta)); end
    for k = 1:numel(ax)
        cla(ax(k)); colorbar(ax(k), 'off'); legend(ax(k), 'off'); title(ax(k), '');
        text(ax(k), 0.5, 0.5, msg, 'Units', 'normalized', 'HorizontalAlignment', 'center');
    end
    return
end
nSel = nnz(sel);
base = staBaseline(R, sta);
[ord, tk, tkLbl, sortLbl] = synOrder(R, G, ftV, sta, sel);
K = R.STA(ord,:) - base(ord);
if numel(ax) >= 4, meanTracePanel(ax(4), R, R.STA - base, sta, sel, selName(sta)); end
cl = kymoClim(K, sta);
imagesc(ax(1), R.t_ms, 1:nSel, K, cl); colormap(ax(1), kymoMap); hold(ax(1), 'on');
xline(ax(1), 0, 'k-'); xline(ax(1), R.score_ms, 'k--');
xlim(ax(1), R.t_ms([1 end])); ylim(ax(1), [0.5 nSel + 0.5]);
if strcmp(sta.synSort, 'x')                              % row ticks: x from the soma (um)
    set(ax(1), 'YTick', tk, 'YTickLabel', tkLbl);
else
    nSig = nnz(R.p(ord) < sta.alpha);
    if nSig > 0 && nSig < nSel, yline(ax(1), nSig + 0.5, 'k-', sprintf('p < %g', sta.alpha), 'LabelHorizontalAlignment', 'left', 'FontSize', 8); end
end
hold(ax(1), 'off');
cb = colorbar(ax(1)); cb.Label.String = gluLabel(sta);
xlabel(ax(1), 'time from spike (ms)'); ylabel(ax(1), sprintf('%s synapse (%s)', selName(sta), sortLbl));
title(ax(1), sprintf('%d spikes; %d %s synapses (spike-tuned %d / %d, %.1f expected)', R.nSpk, nSel, selName(sta), ...
    nnz(R.p < sta.alpha), G.nSyn, sta.alpha*G.nSyn), 'FontSize', 9);
amp = R.amp; amp(~sel) = NaN;                            % unselected footprints are not drawn
if isempty(sta.ampClim)
    a = prctile(amp(sel), [20 90]); if ~(diff(a) > 0), a = [min(amp(sel)) max(amp(sel)) + eps]; end
else
    a = sta.ampClim;
end
footprintMap(ax(2), G, amp, a, sequentialMap, ftV);
title(ax(2), sprintf('amplitude [%g %g] ms (%s synapses)', R.score_ms, selName(sta)), 'FontSize', 9);
cb = colorbar(ax(2)); cb.Label.String = 'z above shuffle'; clim(ax(2), a);
lp = -log10(R.p); pm = max(2, -log10(1/(1 + R.nNull)));
footprintMap(ax(3), G, lp, [0 pm], sequentialMap, ftV, G.xy(sel,:));
title(ax(3), sprintf('-log_{10} p (circles: %s; orange: voltage ROI)', selName(sta)), 'FontSize', 9);
cb = colorbar(ax(3)); cb.Label.String = '-log_{10} p'; clim(ax(3), [0 pm]);
colormap(ax(2), sequentialMap); colormap(ax(3), sequentialMap);
end

function [ord, tk, tkLbl, sortLbl] = synOrder(R, G, ftV, sta, sel)
% Row order of the synapses sel (sta.synSort: 'x' = camera right/basal first, 'p' = lowest p first),
% and for 'x' the row ticks labelled by x distance from the soma (voltage footprint centroid) in um
idx = find(sel);
if strcmp(sta.synSort, 'x')
    [~, o] = sort(G.xy(idx,1), 'descend');
    sortLbl = 'basal \rightarrow apical, x from soma (\mum, + = basal)';
else
    [~, o] = sort(R.p(idx), 'ascend'); sortLbl = 'sorted by p';
end
ord = idx(o);
tk = unique(round(linspace(1, numel(ord), min(numel(ord), 8))));
[yy, xx] = ndgrid(1:size(ftV,1), 1:size(ftV,2)); w = double(ftV);
xSoma = sum(xx(:) .* w(:)) / max(sum(w(:)), eps);
tkLbl = compose('%+.0f', (G.xy(ord(tk), 1) - xSoma) * sta.pixel_um);
end

function base = staBaseline(R, sta)
% Per-synapse baseline of a kymograph: mean over sta.base_ms; sta.base_ms = [] gives 0 (no
% subtraction). When no sample of the kymograph falls in sta.base_ms (window outside sta.kymo_ms),
% the samples before the test window are used instead.
if isempty(sta.base_ms), base = zeros(size(R.STA, 1), 1); return; end
in = R.t_ms >= sta.base_ms(1) & R.t_ms <= sta.base_ms(2);
if ~any(in)
    in = R.t_ms < R.score_ms(1);
    warning('STAglu:baseline', 'sta.base_ms [%g %g] ms is outside the kymograph [%g %g] ms; baseline = %g to %g ms instead', ...
        sta.base_ms, R.t_ms([1 end]), R.t_ms(1), R.score_ms(1));
end
if any(in), base = mean(R.STA(:, in), 2, 'omitnan'); else, base = zeros(size(R.STA, 1), 1); end
end

function cl = kymoClim(K, sta)
% Colour limits of a glutamate kymograph: sta.kymoClim, or the sta.kymoPct percentiles of K when empty
if ~isempty(sta.kymoClim), cl = sta.kymoClim; return; end
pct = [1 99]; if isfield(sta, 'kymoPct') && ~isempty(sta.kymoPct), pct = sta.kymoPct; end
cl = prctile(K(isfinite(K)), pct);
if ~(numel(cl) == 2 && diff(cl) > 0), cl = [-1 1]; end
cl = cl(:)';
end

function s = gluLabel(sta)
% Axis label of the spike-triggered glutamate (with or without baseline subtraction)
if isempty(sta.base_ms), s = 'glutamate z'; else, s = sprintf('glutamate z - baseline [%g %g] ms', sta.base_ms); end
end

function s = selName(sta)
% Name of the synapse set shown in the figures
switch sta.synSelect
    case 'spike', s = 'spike-tuned';
    case 'place', s = 'place-tuned';
    case 'both',  s = 'spike- & place-tuned';
    otherwise,    s = 'all';
end
end

function sel = synSelection(R, G, sta)
% Synapses shown in the figures (sta.synSelect)
spikeT = R.p(:) < sta.alpha; placeT = G.SIp(:) < sta.siAlpha;
switch sta.synSelect
    case 'spike', sel = spikeT;
    case 'place', sel = placeT;
    case 'both',  sel = spikeT & placeT;
    otherwise,    sel = true(G.nSyn, 1);
end
end

function G = addPlaceTuning(G, V, vrFile, sta)
% Spatial information and its p-value for every synapse
S = load(vrFile, 'VR_Result');
[G.SI, G.SIp, G.nEv] = gluSpatialInfo(G, V.t_ax, S.VR_Result, sta);
end

function [SI, SIp, nEv] = gluSpatialInfo(G, tV, VR, sta)
% Place tuning of every synapse, the method of si_pvalue.m / spatialInfo_mat.m (Analysis_20260715 S5):
%  events  : find_spike_bh(zscore(trace), sta.siEventThr) - zscore of the z-scored trace equals
%            zscore(dF/F), so these are the events of S5
%  position: VR_Result (voltage frames) at the voltage frame nearest each glutamate frame; bins of
%            sta.trackLength/sta.placeBins; frames with velocity < sta.velThr excluded
%  SI      : Skaggs information, bits/event = sum_i p_i (r_i/r) log2(r_i/r), p_i occupancy of bin i,
%            r_i event rate in bin i, r the mean rate
%  p       : sta.siShuffle circular shifts of the events (rng 0); p = (1 + #null >= SI) / (1 + nShuffle)
% Default: only running, binned frames are used - they are joined into one running-only timeline,
% r is the mean rate over them, and the events are rolled along that timeline (>= 1 s of running).
% Every shuffle then keeps the synapse's number of running events. SI (bits/event) is strongly
% biased upward for few events (~ (bins-1) / (2 N ln 2)), so when the events are rolled over the
% whole recording (si_pvalue) a shuffle moves events between running and rest; a synapse more active
% at rest than while running then looks place-tuned (null check: 13.5% p < 0.05 instead of 5%).
% sta.siLegacy = true reproduces si_pvalue + spatialInfo_mat exactly (rate over ALL frames, shifts
% over the whole recording; checked: identical SI and p).
ev = find_spike_bh(zscore(G.X, 0, 2), sta.siEventThr(1), sta.siEventThr(2));
[syn, fr] = find(ev); nEv = accumarray(syn, 1, [G.nSyn 1]);
dtV = median(diff(tV));
vIdx = min(max(round((G.t - tV(1)) / dtV) + 1, 1), size(VR, 2));
binTr = ceil(VR(5, vIdx) / (sta.trackLength / sta.placeBins));
binTr(VR(end, vIdx) < sta.velThr) = NaN;
valid = ~isnan(binTr);
if sta.siLegacy
    nB = max(binTr);
    inB = valid & binTr >= 1 & binTr <= nB;
    binCnt = accumarray(binTr(inB)', 1, [nB 1])';
    PR = binCnt / nnz(valid);
    MF = nEv / G.nT;                                     % rate over ALL frames (spatialInfo_mat)
    nT = G.nT; e = fr; bOf = binTr; useE = inB;          % events on the full timeline
else
    nB = sta.placeBins;
    run = find(valid & binTr >= 1 & binTr <= nB);        % running, binned frames
    nT = numel(run);
    tl = zeros(1, G.nT); tl(run) = 1:nT;                 % position of each frame in the running timeline
    e = tl(fr)'; keepE = e > 0;
    syn = syn(keepE); e = e(keepE);
    bOf = binTr(run); useE = true(1, nT);
    binCnt = accumarray(bOf', 1, [nB 1])';
    PR = binCnt / nT;
    MF = accumarray(syn, 1, [G.nSyn 1]) / nT;            % running events / running frames (shift invariant)
end
SI = siOfEvents(e);
minShift = min(max(1, round(1/G.dt)), max(1, floor(nT/2)));
rng(0); shifts = randi([minShift, max(minShift+1, nT-minShift)], 1, sta.siShuffle);
nGE = zeros(G.nSyn, 1);
for k = 1:sta.siShuffle
    nGE = nGE + (siOfEvents(mod(e - 1 + shifts(k), nT) + 1) >= SI);
end
SIp = (1 + nGE) / (1 + sta.siShuffle);
SIp(isnan(SI)) = NaN;
    function si = siOfEvents(frames)
        bb = bOf(frames); kk = useE(frames);
        binSum = accumarray([syn(kk) bb(kk)'], 1, [G.nSyn nB]);
        FR = binSum ./ binCnt;
        si = sum(PR .* FR .* log2(FR ./ MF), 2, 'omitnan') ./ MF;
    end
end

function meanTracePanel(ax, R, K, sta, sel, selLbl)
% Mean (+- SEM) baseline-subtracted spike-triggered glutamate of the selected synapses and of all
% (when all are selected: the p < alpha synapses and all)
cla(ax); hold(ax, 'on');
if all(sel), grp = {R.p(:) < sta.alpha, sel}; lbl1 = sprintf('p < %g', sta.alpha);
else,        grp = {sel, true(size(sel))};    lbl1 = selLbl; end
col = {[0.922 0.408 0.204], [0.165 0.471 0.839]};
lbl = {sprintf('%s (%d)', lbl1, nnz(grp{1})), sprintf('all synapses (%d)', numel(sel))};
h = gobjects(0); L = {};
for g = 1:2
    X = double(K(grp{g}, :));
    if isempty(X), continue; end
    m = mean(X, 1, 'omitnan'); se = std(X, 0, 1, 'omitnan') / sqrt(size(X,1));
    fill(ax, [R.t_ms fliplr(R.t_ms)], [m + se, fliplr(m - se)], col{g}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    h(end+1) = plot(ax, R.t_ms, m, 'Color', col{g}, 'LineWidth', 2); L{end+1} = lbl{g}; %#ok<AGROW>
end
xline(ax, 0, 'k-'); xline(ax, R.score_ms, 'k--'); yline(ax, 0, ':', 'Color', [0.5 0.5 0.5]);
hold(ax, 'off'); xlim(ax, R.t_ms([1 end])); box(ax, 'off');
xlabel(ax, 'time from spike (ms)'); ylabel(ax, gluLabel(sta));
legend(ax, h, L, 'Location', 'northwest', 'Box', 'off', 'FontSize', 8);
title(ax, 'mean spike-triggered glutamate (\pm SEM)', 'FontSize', 9);
end

function footprintMap(ax, G, v, lim, cmap, ftV, mark)
% Glutamate mean image with every footprint drawn with its own weights (each scaled to peak 1) in the
% colour of its value v; where footprints overlap the colours are weight-averaged and the opacity is
% the summed weight (capped at 1). The voltage neuron's footprint is outlined orange; mark ([x y])
% circled in [0 0.6 1].
img = mat2gray(G.img, double(prctile(G.img(:), [1 99.7])));
bg = repmat(img * 0.7, [1 1 3]);
v = v(:); ok = isfinite(v);
ci = round(1 + (size(cmap,1)-1) * min(max((v(ok) - lim(1)) / diff(lim), 0), 1));
S = G.S(:, ok);
S = S * spdiags(1 ./ max(full(max(S, [], 1)), eps)', 0, size(S,2), size(S,2));   % each footprint peak 1
A = full(S * ones(size(S,2), 1));                        % summed weight per pixel
Cw = full(S * cmap(ci,:)) ./ max(A, eps);                % weight-averaged colour per pixel
al = min(A, 1);
rgb = bg;
for c = 1:3
    rgb(:,:,c) = reshape((1 - al) .* reshape(bg(:,:,c), [], 1) + al .* Cw(:,c), size(img));
end
image(ax, rgb); axis(ax, 'image', 'off'); hold(ax, 'on');
if ~isempty(ftV) && any(ftV(:) > 0)
    contour(ax, ftV, max(ftV(:)) * [0.5 0.5], 'EdgeColor', [1 0.55 0], 'LineWidth', 1.25);   % voltage ROI: orange
end
if nargin > 6 && ~isempty(mark)
    plot(ax, mark(:,1), mark(:,2), 'o', 'Color', [0 0.6 1], 'MarkerSize', 4, 'LineWidth', 0.75);   % glutamate markers
end
hold(ax, 'off');
end

function c = kymoMap
% Kymograph colormap: blue - white - red
c = gen_colormap([0 0.5 1; 1 1 1; 1 0 0], 256);
end

function c = sequentialMap
% dark gray -> orange -> light yellow
k = [0.2 0.2 0.2; 0.922 0.408 0.204; 1 0.95 0.6];
c = interp1([0 0.6 1], k, linspace(0, 1, 256));
end

function [keep, preISI, postISI] = selectSpikes(V, n, type, preMin_ms, postMin_ms)
% Spikes of ROI n of one type ('all', 'simple' = SpClass(1), 'complex' = SpClass(2) onsets) whose
% pre-ISI (since the previous spike) >= preMin_ms and post-ISI (to the next spike; from the end of
% the burst for a complex-spike onset) >= postMin_ms. Returns voltage frame indices.
spAll = find(V.spike(n,:) > 0);
switch type
    case 'simple',  cand = find(V.SpClass(1,:,n) > 0);
    case 'complex', cand = find(V.SpClass(2,:,n) > 0);
    otherwise,      cand = spAll;
end
t = V.t_ax(:)' * 1000;                                   % ms
lastFr = cand;                                           % last frame of the event (burst end for CS)
if strcmp(type, 'complex') && isfield(V, 'CStrace')
    L = bwlabel(V.CStrace(n,:) > 0);
    e = accumarray(L(L > 0)', find(L > 0)', [], @max);
    hasL = L(cand) > 0;
    lastFr(hasL) = max(cand(hasL), e(L(cand(hasL)))');
end
ip = discretize(cand, [-Inf spAll Inf]) - 1;             % index of the last spike at or before cand
prevIdx = ip - (spAll(max(ip,1)) == cand);               % previous spike strictly before
preISI = inf(size(cand)); ok = prevIdx >= 1;
preISI(ok) = t(cand(ok)) - t(spAll(prevIdx(ok)));
in = discretize(lastFr, [-Inf spAll Inf]);               % first spike strictly after the event end
postISI = inf(size(cand)); ok = in <= numel(spAll);
postISI(ok) = t(spAll(in(ok))) - t(lastFr(ok));
keep = cand(preISI >= preMin_ms & postISI >= postMin_ms);
end

function interactiveSpikeTrigGlu(V, n, G, sta, ttl)
% Docked figure: spike type, pre-/post-ISI minimum and glutamate window; Update reruns the selection,
% the kymograph, the permutation test and the maps
fig = figure('Name', ['Spike-triggered glutamate ' ttl], 'NumberTitle', 'off', 'WindowStyle', 'docked');   % UI components in a dockable figure
try, theme(fig, 'light'); catch, end
gl = uigridlayout(fig, [1 2], 'ColumnWidth', {230, '1x'});
cp = uigridlayout(gl, [14 1], 'RowHeight', [repmat({28}, 1, 13), {'1x'}]);
uilabel(cp, 'Text', ttl, 'FontWeight', 'bold');
uilabel(cp, 'Text', 'Spike type');
dd = uidropdown(cp, 'Items', {'All spikes', 'Simple spikes', 'Complex-spike onsets'}, 'ItemsData', {'all', 'simple', 'complex'});
uilabel(cp, 'Text', 'Glutamate window start / end (ms)');
wg = uigridlayout(cp, [1 2], 'Padding', [0 0 0 0]);
w1 = uieditfield(wg, 'numeric', 'Value', sta.score_ms(1)); w2 = uieditfield(wg, 'numeric', 'Value', sta.score_ms(2));
uilabel(cp, 'Text', 'Pre-ISI >= (ms)');  pre  = uieditfield(cp, 'numeric', 'Value', 0, 'Limits', [0 Inf]);
uilabel(cp, 'Text', 'Post-ISI >= (ms)'); post = uieditfield(cp, 'numeric', 'Value', 0, 'Limits', [0 Inf]);
uilabel(cp, 'Text', 'Synapses shown');
ddS = uidropdown(cp, 'Items', {sprintf('Spike-tuned (p < %g)', sta.alpha), sprintf('Place-tuned (SI p < %g)', sta.siAlpha), ...
    'Spike- & place-tuned', 'All'}, 'ItemsData', {'spike', 'place', 'both', 'all'}, 'Value', sta.synSelect);
ddO = uidropdown(cp, 'Items', {'Sort: x, right (basal) to left', 'Sort: p'}, 'ItemsData', {'x', 'p'}, 'Value', sta.synSort);
btn = uibutton(cp, 'Text', 'Update', 'FontWeight', 'bold');
st = uilabel(cp, 'Text', '', 'WordWrap', 'on');
pg = uigridlayout(gl, [3 2], 'RowHeight', {'1x', '1x', '1.1x'});
ax = [uiaxes(pg) uiaxes(pg) uiaxes(pg) uiaxes(pg)];
axV = uiaxes(pg);
ax(1).Layout.Row = [1 2]; ax(1).Layout.Column = 1;          % kymograph
ax(4).Layout.Row = 1;     ax(4).Layout.Column = 2;          % mean glutamate traces
axV.Layout.Row = 2;       axV.Layout.Column = 2;            % spike-triggered voltage
ax(2).Layout.Row = 3;     ax(2).Layout.Column = 1;          % amplitude map
ax(3).Layout.Row = 3;     ax(3).Layout.Column = 2;          % p-value map
trV = voltZ(V, n);
dtV = median(diff(V.t_ax)) * 1000;
btn.ButtonPushedFcn = @(~,~) update();
update();
    function update()
        btn.Enable = 'off'; st.Text = 'computing ...'; drawnow;
        keep = selectSpikes(V, n, dd.Value, pre.Value, post.Value);
        win = sort([w1.Value w2.Value]);
        R = spikeTrigGlu(G, V.t_ax(keep), win, sta);
        s2 = sta; s2.synSelect = ddS.Value; s2.synSort = ddO.Value;
        voltTracePanel(axV, trV, keep(R.used), dtV, sta, win);
        spikeTrigGluPanels(R, G, V.ftprnt_glu(:,:,n), s2, ax);
        st.Text = sprintf('%d spikes used; spike-tuned %d / %d synapses (%.1f expected by chance), %d of them place-tuned; place-tuned %d.', ...
            R.nSpk, nnz(R.p < sta.alpha), G.nSyn, sta.alpha*G.nSyn, nnz(R.p < sta.alpha & G.SIp < sta.siAlpha), nnz(G.SIp < sta.siAlpha));
        assignin('base', 'STAgluInteractive', R);        % last result in the workspace
        btn.Enable = 'on';
    end
end

function S = loadRetry(fn, varargin)
% load() for -v7.3 files on the network share; on an SMB/HDF5 error, load a local copy instead
try
    S = load(fn, varargin{:}); return
catch ME
    fprintf('  network load failed (%s) -> loading a local copy\n', ME.message);
end
tmp = [tempname '.mat'];
copyfile(fn, tmp);
c = onCleanup(@() delete(tmp));
S = load(tmp, varargin{:});
end

function interactiveEventViewer(V, n, G, sta, ttl)
% Docked figure stepping through the events of ROI n (left/right arrows, Prev/Next, or the event
% number): voltage around the event, and glutamate kymographs of the tuned / other synapses
% (rows basal = camera right at the top -> apical, labelled in um from the soma)
fig = figure('Name', ['Event viewer ' ttl], 'NumberTitle', 'off', 'WindowStyle', 'docked');
try, theme(fig, 'light'); catch, end
gl = uigridlayout(fig, [1 2], 'ColumnWidth', {230, '1x'});
cp = uigridlayout(gl, [13 1], 'RowHeight', [repmat({28}, 1, 12), {'1x'}]);
uilabel(cp, 'Text', ttl, 'FontWeight', 'bold');
uilabel(cp, 'Text', 'Event type');
dd = uidropdown(cp, 'Items', {'Complex-spike onsets', 'Simple spikes', 'All spikes'}, 'ItemsData', {'complex', 'simple', 'all'});
uilabel(cp, 'Text', 'Isolation: pre-ISI / post-ISI >= (ms)');
ig = uigridlayout(cp, [1 2], 'Padding', [0 0 0 0]);
preE = uieditfield(ig, 'numeric', 'Value', 0, 'Limits', [0 Inf]); postE = uieditfield(ig, 'numeric', 'Value', 0, 'Limits', [0 Inf]);
uilabel(cp, 'Text', 'Event (left/right arrow keys)');
ng = uigridlayout(cp, [1 3], 'Padding', [0 0 0 0], 'ColumnWidth', {50, '1x', 50});
bPrev = uibutton(ng, 'Text', '<'); ev = uieditfield(ng, 'numeric', 'Value', 1, 'RoundFractionalValues', 'on'); bNext = uibutton(ng, 'Text', '>');
cbAvg = uicheckbox(cp, 'Text', 'Average of all events');
uilabel(cp, 'Text', 'Tuned synapses');
ddT = uidropdown(cp, 'Items', {sprintf('Place-tuned (SI p < %g)', sta.siAlpha), sprintf('Spike-triggered (p < %g)', sta.alpha), 'Both'}, ...
    'ItemsData', {'place', 'spike', 'both'});
uilabel(cp, 'Text', 'Window (ms) / colour percentiles');
pg2 = uigridlayout(cp, [1 4], 'Padding', [0 0 0 0], 'ColumnSpacing', 4);
wA = uieditfield(pg2, 'numeric', 'Value', -500); wB = uieditfield(pg2, 'numeric', 'Value', 500);
pk = [1 99]; if isfield(sta, 'kymoPct') && ~isempty(sta.kymoPct), pk = sta.kymoPct; end
pLo = uieditfield(pg2, 'numeric', 'Value', pk(1), 'Limits', [0 100]); pHi = uieditfield(pg2, 'numeric', 'Value', pk(2), 'Limits', [0 100]);
st = uilabel(cp, 'Text', '', 'WordWrap', 'on', 'VerticalAlignment', 'top');
rg = uigridlayout(gl, [3 1]);
axV = uiaxes(rg); axT = uiaxes(rg); axU = uiaxes(rg);
trV = voltZ(V, n); tV = V.t_ax(:)'; dtV = median(diff(tV));
spAll = find(V.spike(n,:) > 0);
if isfield(V, 'CStrace'), Lcs = bwlabel(V.CStrace(n,:) > 0); else, Lcs = zeros(1, size(V.spike,2)); end
cache = struct();                                        % spike-triggered results per event type
events = []; type = ''; nAll = 0;
dd.ValueChangedFcn = @(~,~) newType();
preE.ValueChangedFcn = @(~,~) newType(); postE.ValueChangedFcn = @(~,~) newType();   % new isolation: reselect
bPrev.ButtonPushedFcn = @(~,~) step(-1); bNext.ButtonPushedFcn = @(~,~) step(1);
ev.ValueChangedFcn = @(~,~) draw();
for hc = [cbAvg ddT wA wB pLo pHi], hc.ValueChangedFcn = @(~,~) draw(); end
ftV = V.ftprnt_glu(:,:,n); [yy, xx] = ndgrid(1:size(ftV,1), 1:size(ftV,2));
xSoma = sum(xx(:) .* double(ftV(:))) / max(sum(double(ftV(:))), eps);   % soma: voltage footprint centroid
linkaxes([axV axT axU], 'x');
fig.WindowKeyPressFcn = @keyPress;
newType();
    function keyPress(~, e)
        switch e.Key
            case 'rightarrow', step(1);
            case 'leftarrow',  step(-1);
        end
    end
    function newType()
        type = dd.Value;
        nAll = numel(selectSpikes(V, n, type, 0, 0));
        events = selectSpikes(V, n, type, preE.Value, postE.Value);   % isolated: pre-/post-ISI limits
        events = events(tV(events) > G.t(1) & tV(events) < G.t(end));   % inside the glutamate recording
        ev.Limits = [1 max(1, numel(events))]; ev.Value = 1;
        draw();
    end
    function step(d)
        if isempty(events), return; end
        ev.Value = min(max(ev.Value + d, 1), numel(events));
        draw();
    end
    function R = stRes()
        % Section 3 spike-triggered result of all events of this type (cached)
        win = sort([wA.Value wB.Value]); iso = [preE.Value postE.Value];
        if ~isfield(cache, type) || ~isequal(cache.(type).win, win) || ~isequal(cache.(type).iso, iso)
            st.Text = 'computing the spike-triggered average ...'; drawnow;   % new type, window or isolation
            s2 = sta; s2.kymo_ms = [min(win(1), -1) max(win(2), 1)];
            cache.(type) = spikeTrigGlu(G, tV(events), sta.score_ms, s2);
            cache.(type).win = win; cache.(type).iso = iso;
        end
        R = cache.(type);
    end
    function draw()
        if isempty(events)
            cla(axV); cla(axT); cla(axU); st.Text = sprintf('no events (%d before the ISI limits)', nAll); return
        end
        win = sort([wA.Value wB.Value]);
        k = ev.Value; e0 = events(k); t0 = tV(e0);
        % tuned set
        switch ddT.Value
            case 'place', tuned = G.SIp < sta.siAlpha;
            case 'spike', R = stRes(); tuned = R.p < sta.alpha;
            otherwise,    R = stRes(); tuned = G.SIp < sta.siAlpha & R.p < sta.alpha;
        end
        % voltage
        o = round(win(1)/1000/dtV):round(win(2)/1000/dtV);
        cla(axV); hold(axV, 'on');
        if cbAvg.Value
            R = stRes();
            kk = events(events + o(1) >= 1 & events + o(end) <= numel(trV));
            plot(axV, o*dtV*1000, mean(trV(kk(:) + o), 1), 'k', 'LineWidth', 1);
            title(axV, sprintf('average of %d %s events', numel(kk), dd.Items{strcmp(dd.ItemsData, type)}), 'FontSize', 10);
        else
            ii = e0 + o; ok = ii >= 1 & ii <= numel(trV);
            plot(axV, o(ok)*dtV*1000, trV(ii(ok)), 'k', 'LineWidth', 0.75);
            sp = spAll(spAll >= ii(find(ok,1)) & spAll <= ii(find(ok,1,'last')));
            plot(axV, (sp - e0)*dtV*1000, trV(sp), '.', 'Color', [0.85 0.1 0.1], 'MarkerSize', 12);
            nSp = 1; if strcmp(type, 'complex') && Lcs(e0) > 0, nSp = nnz(ismember(spAll, find(Lcs == Lcs(e0)))); end
            [~, gk] = min(abs(G.t - t0));
            title(axV, sprintf('Event %d / %d: ROI %d %s (nSpike = %d)   Vframe %d / Gframe %d', k, numel(events), n, ...
                dd.Items{strcmp(dd.ItemsData, type)}, nSp, e0, gk), 'FontSize', 10);
        end
        xline(axV, 0, 'Color', [0.85 0.1 0.1]); hold(axV, 'off'); xlim(axV, win); ylabel(axV, 'voltage (robust z)'); box(axV, 'off');
        % glutamate kymographs: tuned / other synapses, rows basal (camera right) -> apical
        if cbAvg.Value
            tg = R.t_ms; Y = R.STA;                       % nSyn x time, interpolated at the voltage rate
        else
            fr = find(G.t >= t0 + win(1)/1000 & G.t <= t0 + win(2)/1000);
            tg = (G.t(fr) - t0) * 1000; Y = G.X(:, fr);   % glutamate frames of this event
        end
        grp = {find(tuned), find(~tuned)}; axs = [axT axU];
        nm = {sprintf('tuned (%s)', ddT.Items{strcmp(ddT.ItemsData, ddT.Value)}), 'other synapses'};
        s2 = sta; s2.kymoClim = []; s2.kymoPct = sort([pLo.Value pHi.Value]);
        for g = 1:2
            ax = axs(g); cla(ax); colorbar(ax, 'off');
            idx = grp{g};
            if isempty(idx) || isempty(tg), title(ax, sprintf('%s: none', nm{g}), 'FontSize', 10); continue; end
            [~, o] = sort(G.xy(idx, 1), 'descend'); idx = idx(o);
            K = double(Y(idx, :));
            imagesc(ax, tg, 1:numel(idx), K, kymoClim(K, s2)); colormap(ax, kymoMap); hold(ax, 'on');
            xline(ax, 0, 'k-', 'LineWidth', 1); hold(ax, 'off');
            tk = unique(round(linspace(1, numel(idx), min(numel(idx), 7))));
            set(ax, 'YDir', 'reverse', 'YTick', tk, 'YTickLabel', compose('%+.0f', (G.xy(idx(tk), 1) - xSoma) * sta.pixel_um));
            ylim(ax, [0.5 numel(idx) + 0.5]);
            ylabel(ax, {'basal \rightarrow apical', '\mum from soma'});
            cb = colorbar(ax); cb.Label.String = 'glutamate z';
            title(ax, sprintf('%s (%d synapses)', nm{g}, numel(idx)), 'FontSize', 10);
        end
        xlabel(axU, 'time from event (ms)');
        st.Text = sprintf(['%d tuned / %d synapses; %d of %d events pass pre-ISI >= %g ms, post-ISI >= %g ms ' ...
            '(complex spikes: post-ISI from the burst end). Arrow keys: previous / next event.'], ...
            nnz(tuned), G.nSyn, numel(events), nAll, preE.Value, postE.Value);
    end
end

function C = classifySpikes(V, n, silent_ms, burstISI_ms)
% Voltage frames of i) isolated spikes, ii) first spikes of bursts (not complex spikes), iii) first
% spikes of complex spikes, for silent window silent_ms (see Section 3c)
spAll = find(V.spike(n,:) > 0);
t = V.t_ax(spAll) * 1000; t = t(:)';
prevISI = [Inf diff(t)]; nextISI = [diff(t) Inf];
inCS = false(size(spAll));
if isfield(V, 'CStrace'), inCS = V.CStrace(n, spAll) > 0; end
C.isolated   = spAll(prevISI >= silent_ms & nextISI >= silent_ms);
C.burstFirst = spAll(prevISI >= silent_ms & nextISI <= burstISI_ms & ~inCS);
C.csFirst    = selectSpikes(V, n, 'complex', silent_ms, 0);
end

function fig = spikeClassKymographs(V, n, G, sta, figNo, ttl)
% Voltage average (top) and glutamate kymograph of the spike-tuned synapses for the three spike
% classes of classifySpikes; rows basal (camera right, top) to apical (left, bottom), in um from the soma
C = classifySpikes(V, n, sta.silent_ms, sta.burstISI_ms);
spAll = find(V.spike(n,:) > 0);
Rall = spikeTrigGlu(G, V.t_ax(spAll), sta.score_ms, sta);     % spike-tuned synapses: all spikes
tuned = find(Rall.p < sta.alpha);
[~, o] = sort(G.xy(tuned, 1), 'descend'); rows = tuned(o);   % camera right (basal) first
ft = double(V.ftprnt_glu(:,:,n)); [yy, xx] = ndgrid(1:size(ft,1), 1:size(ft,2));
xSoma = sum(xx(:) .* ft(:)) / sum(ft(:));
dx_um = (G.xy(rows, 1) - xSoma) * sta.pixel_um;            % + = basal (right of the soma)
cls = {C.isolated, C.burstFirst, C.csFirst};
nm = {sprintf('i) isolated (\\pm%g ms)', sta.silent_ms), sprintf('ii) burst first (ISI \\leq %g ms)', sta.burstISI_ms), ...
      'iii) complex-spike first'};
R = cell(1, 3); K = cell(1, 3);
for c = 1:3
    R{c} = spikeTrigGlu(G, V.t_ax(cls{c}), sta.score_ms, sta);
    base = staBaseline(R{c}, sta);
    K{c} = R{c}.STA(rows, :) - base(rows);
end
allK = cell2mat(cellfun(@(k) k(:), K, 'UniformOutput', false)');
cl = kymoClim(allK, sta);
fig = figure(figNo); clf; set(fig, 'Name', [ttl ' spike classes']);
try, theme(fig, 'light'); catch, set(fig, 'Color', 'w'); end
tl = tiledlayout(fig, 4, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
trV = voltZ(V, n); dtV = median(diff(V.t_ax)) * 1000;
tk = unique(round(linspace(1, numel(rows), min(numel(rows), 9))));
axAll = gobjects(0);
for c = 1:3
    ax = nexttile(tl, c);                                % voltage average of the class
    voltTracePanel(ax, trV, cls{c}(R{c}.used), dtV, sta, sta.score_ms); axAll(end+1) = ax; %#ok<AGROW>
    title(ax, sprintf('%s, n = %d', nm{c}, R{c}.nSpk), 'FontSize', 9); xlabel(ax, '');
    ax = nexttile(tl, 3 + c, [3 1]); axAll(end+1) = ax; %#ok<AGROW>
    if R{c}.nSpk == 0 || isempty(rows)
        text(ax, 0.5, 0.5, 'no spikes / no spike-tuned synapses', 'Units', 'normalized', 'HorizontalAlignment', 'center'); axis(ax, 'off'); continue
    end
    imagesc(ax, R{c}.t_ms, 1:numel(rows), K{c}, cl); colormap(ax, kymoMap); hold(ax, 'on');
    xline(ax, 0, 'k-'); xline(ax, sta.score_ms, 'k--'); hold(ax, 'off');
    set(ax, 'YTick', tk, 'YTickLabel', compose('%+.0f', dx_um(tk)));
    xlabel(ax, 'time from spike (ms)');
    if c == 1, ylabel(ax, {'spike-tuned synapses (basal \rightarrow apical)', 'x from soma (\mum, + = basal)'}); end
    if c == 3, cb = colorbar(ax); cb.Label.String = gluLabel(sta); end
end
linkaxes(axAll, 'x'); xlim(axAll(1), sta.kymo_ms);   % all panels share the time axis (zoom/pan together)
title(tl, sprintf('%s: %d spike-tuned synapses (all spikes, p < %g); silent window %g ms', ttl, numel(rows), ...
    sta.alpha, sta.silent_ms), 'FontSize', 10);
drawnow;
end

function fig = representativeFigure(V, G, VR, rep, top, synLbl, figNo, ttl)
% Section 3a0 figure: voltage image + ROI of neuron rep.n, glutamate image + footprints of the
% synapses top in the given order (described by synLbl; hsv by row, circled), voltage trace of rep.n with
% spikes and complex spikes, glutamate kymograph of top (rows in that order, row colours marked at the
% left as the footprints), VR position
n = rep.n; cOr = [1 0.55 0];
top = top(:);
cT = hsv(numel(top));                                    % row colours (footprints, circles, row bullets)
tV = V.t_ax(1:size(V.traces, 2)); tV = tV(:)';
if isempty(rep.win_s), win = [max(tV(1), G.t(1)) min(tV(end), G.t(end))]; else, win = rep.win_s; end

fig = figure(figNo); clf; set(fig, 'Name', ttl);
try, theme(fig, 'light'); catch, set(fig, 'Color', 'w'); end
tl = tiledlayout(fig, 8, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% voltage reference image in the glutamate frame (cropped to the voltage camera FOV) + ROI of neuron n
ax = nexttile(tl, 1, [2 1]);
R = V.ref_im_glu; v = R(R > 0);
[rr, cc] = find(R > 0); r0 = min(rr); c0 = min(cc);
imshow(R(r0:max(rr), c0:max(cc)), prctile(v, [1 99.7]), 'Parent', ax); hold(ax, 'on');
P = V.ROIpoly_glu{n};
plot(ax, P([1:end 1],1) - c0 + 1, P([1:end 1],2) - r0 + 1, '-', 'Color', cOr, 'LineWidth', 1.5);
hold(ax, 'off'); title(ax, sprintf('voltage (Voltron2-ST) in the glutamate frame; ROI %d', n), 'FontSize', 9);
% glutamate mean image + footprints of the most active synapses
ax = nexttile(tl, 2, [2 1]);
imshow(G.img, prctile(G.img(:), [1 99.7]), 'Parent', ax); hold(ax, 'on');
[H, W] = size(G.img);
for k = 1:numel(top)
    m = reshape(full(G.S(:, top(k))), H, W);
    b = bwboundaries(m > 0.5 * max(m(:)));
    for j = 1:numel(b), plot(ax, b{j}(:,2), b{j}(:,1), '-', 'Color', cT(k,:), 'LineWidth', 1); end
    plot(ax, G.xy(top(k),1), G.xy(top(k),2), 'o', 'Color', cT(k,:), 'MarkerSize', rep.markerSize, 'LineWidth', 2);   % circle: visible
end
plot(ax, P([1:end 1],1), P([1:end 1],2), '-', 'Color', cOr, 'LineWidth', 1);
hold(ax, 'off'); title(ax, sprintf('glutamate (iGluSnFR) mean; %s (colour: row)', synLbl), 'FontSize', 9);
% voltage trace of neuron n: spikes and complex spikes
axT = gobjects(0);
ax = nexttile(tl, 5, [1 2]); axT(end+1) = ax; hold(ax, 'on');
iv = find(tV >= win(1) & tV <= win(2));
z = voltZ(V, n);
yl = [min(z(iv)) max(z(iv))] + [-0.5 0.5];
hc = gobjects(0);
if isfield(V, 'CStrace')                                 % complex spikes: shaded
    cs = V.CStrace(n, iv) > 0; d = diff([0 cs 0]); on = find(d == 1); off = find(d == -1) - 1;
    for k = 1:numel(on)
        hc = patch(ax, tV(iv([on(k) off(k) off(k) on(k)])), yl([1 1 2 2]), cOr, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    end
end
plot(ax, tV(iv), z(iv), 'k-', 'LineWidth', 0.5);
sp = iv(V.spike(n, iv) > 0);
hs = plot(ax, tV(sp), z(sp), '.', 'Color', [0.85 0.1 0.1], 'MarkerSize', 8);
hold(ax, 'off'); ylim(ax, yl); ylabel(ax, 'voltage (z)'); box(ax, 'off');
L = {sprintf('spikes (%d)', numel(sp))}; h = hs;
if ~isempty(hc), h = [hs hc]; L{end+1} = sprintf('complex spikes (%d)', numel(on)); end
legend(ax, h, L, 'Location', 'northeast', 'Box', 'off', 'FontSize', 8);
title(ax, sprintf('voltage, ROI %d (robust z)', n), 'FontSize', 9);
% glutamate kymograph (rows in the given order; tick marks coloured as the footprints)
ax = nexttile(tl, 7, [4 2]); axT(end+1) = ax;
ig = find(G.t >= win(1) & G.t <= win(2));
Kz = double(G.X(top, ig));
% G.X is already the robust z of each trace: (x - median) / (1.4826 MAD), over the whole recording
cl = rep.clim; if isempty(cl), c = prctile(Kz(:), 99.5); if ~(c > 0), c = 1; end, cl = [0 c]; end
imagesc(ax, G.t(ig([1 end])), [1 numel(top)], Kz, cl); colormap(ax, kymoMap);
cb = colorbar(ax, 'Location', 'east'); cb.Label.String = 'glutamate (robust z, MAD)';   % inside: layout and time axes unchanged
set(ax, 'YTick', 1:numel(top), 'TickLabelInterpreter', 'tex', 'TickLength', [0 0], ...
    'YTickLabel', arrayfun(@(k) sprintf('\\fontsize{14}\\color[rgb]{%.3f,%.3f,%.3f}\\bullet', cT(k,:)), 1:numel(top), 'UniformOutput', false));
ylabel(ax, {sprintf('%d spike-tuned synapses', numel(top)), 'basal (top) \rightarrow apical'});
box(ax, 'off');
% VR position + reward
ax = nexttile(tl, 15, [1 2]); axT(end+1) = ax; hold(ax, 'on');
nVR = min(size(VR, 2), numel(tV)); tvr = tV(1:nVR); ir = find(tvr >= win(1) & tvr <= win(2));
plot(ax, tvr(ir), VR(5, ir), 'k-', 'LineWidth', 0.75);
rw = ir(VR(11, ir) > 0 & [true, VR(11, ir(1:end-1)) == 0]);   % reward onsets
plot(ax, tvr(rw), VR(5, rw), 'v', 'Color', [0.165 0.471 0.839], 'MarkerFaceColor', [0.165 0.471 0.839], 'MarkerSize', 6);
hold(ax, 'off'); ylabel(ax, 'position'); xlabel(ax, 'time (s)'); box(ax, 'off');
title(ax, 'VR track position (triangles: reward)', 'FontSize', 9);
linkaxes(axT, 'x'); xlim(axT(1), win);
for a = axT(1:end-1), set(a, 'XTickLabel', []); end
title(tl, ttl, 'Interpreter', 'none', 'FontSize', 10);
drawnow;
end

function V = V_forMovie(fp)
% Volt_Result fields needed by Section 3d
S = load(fullfile(fp, 'Volt_Result.mat'), 'Result'); R = S.Result;
V = struct('spike', R.spike, 'SpClass', R.SpClass, 't_ax', R.t_ax, 'ROIpoly', {R.ROIpoly}, 'ROIpoly_glu', {R.ROIpoly_glu});
V.CStrace = zeros(size(R.spike)); if isfield(R, 'CStrace'), V.CStrace = R.CStrace; end
end

function SM = trigClassMovies(fp, V, n, sta, mov)
% Spike-triggered average voltage and glutamate movies of three spike classes from the snippets
od = fullfile(fp, 'TriggerSnippets');
X = load(fullfile(od, 'Trig_index.mat'), 'index', 'common');
ii = find([X.index.roi] == n, 1);
if isempty(ii), error('ROI %d has no snippets in %s (only good neurons are saved)', n, od); end
% saved events: voltage frame -> (bin, event)
evF = []; evBin = []; evIdx = []; evShift = []; evT = [];
for b = 1:numel(X.index(ii).bins)
    S = load(fullfile(od, X.index(ii).bins{b}), 'info');
    k = S.info.nEvents;
    evF = [evF, S.info.events.voltFrame]; evBin = [evBin, b*ones(1,k)]; evIdx = [evIdx, 1:k]; %#ok<AGROW>
    evShift = [evShift, S.info.events.tShift_s]; evT = [evT, S.info.events.t_s]; %#ok<AGROW>
end
% classes from the current spikes (as Section 3c)
C = classifySpikes(V, n, sta.silent_ms, sta.burstISI_ms);
cls = {C.isolated, C.burstFirst, find(V.SpClass(2,:,n) > 0)};
SM.names = {sprintf('isolated (\\pm%g ms)', sta.silent_ms), sprintf('burst first (ISI \\leq %g ms)', sta.burstISI_ms), ...
    'complex-spike first'};
gI = X.common.glu; vI = X.common.volt;
lagFit = [0 0]; if isfield(sta, 'lagCorrect') && sta.lagCorrect, lagFit = gluLagOf(fp, sta); end
SM.tV = vI.frameTimes_s; SM.tG = gI.frameTimes_s;        % voltage frames / glutamate grid (s from the spike)
dtG = 1 / gI.frameRate;
for c = 1:3
    [ok, loc] = ismember(cls{c}, evF);
    SM.nClass(c) = numel(cls{c}); SM.nEvents(c) = nnz(ok);
    loc = loc(ok);
    fprintf('  %s: %d of %d events in the snippets\n', SM.names{c}, numel(loc), numel(cls{c}));
    if numel(loc) > mov.maxEvents, loc = loc(round(linspace(1, numel(loc), mov.maxEvents))); end
    SM.nEvents(c) = numel(loc);
    sumV = zeros([vI.frameSize vI.nFramesPerEvent], 'single'); sumG = zeros([gI.frameSize gI.nFramesPerEvent], 'single');
    cntG = zeros(1, gI.nFramesPerEvent);
    tc = tic;
    for e = 1:numel(loc)
        mf = fullfile(od, X.index(ii).bins{evBin(loc(e))});
        v = readTrigSnippets(mf, 'volt', evIdx(loc(e)));
        sumV = sumV + v;
        g = readTrigSnippets(mf, 'glu', evIdx(loc(e)));
        ts = SM.tG - evShift(loc(e)) + polyval(lagFit, evT(loc(e)))/1000 + sta.frameMid * dtG/2;   % frame times from the spike (+ lag if sta.lagCorrect)
        D = reshape(g, [], numel(ts))';
        Di = interp1(ts(:), D, SM.tG(:), 'linear');         % onto the common grid; NaN outside
        okF = all(isfinite(Di) | isnan(D(1,:)), 2)';         % grid frames covered
        Di(~isfinite(Di)) = 0;
        sumG = sumG + reshape(Di', size(sumG)); cntG = cntG + okF;
        if mod(e, 10) == 0, fprintf('    %d / %d events (%.0f s)\n', e, numel(loc), toc(tc)); end
    end
    mF = sumV / max(SM.nEvents(c), 1);
    SM.volt{c} = -(mF - median(mF, 3));                      % -dF vs the pixel's temporal median (Voltron dims)
    mF = sumG ./ reshape(max(cntG, 1), 1, 1, []); mF(:,:,cntG == 0) = NaN;
    SM.glu{c} = mF - median(mF, 3, 'omitnan'); SM.glu{c}(isnan(SM.glu{c})) = 0;
    if mov.sigma > 0, SM.volt{c} = imgaussfilt(SM.volt{c}, mov.sigma); SM.glu{c} = imgaussfilt(SM.glu{c}, mov.sigma); end
end
SM.roiV = V.ROIpoly{n}; SM.roiG = V.ROIpoly_glu{n};
mV = poly2mask(SM.roiV(:,1), SM.roiV(:,2), vI.frameSize(1), vI.frameSize(2));
for c = 1:3                                              % voltage trace: movie mean in the ROI
    X = reshape(SM.volt{c}, [], numel(SM.tV)); SM.voltTrace{c} = mean(X(mV(:), :), 1);
end
end

function fig = trigMovieViewer(SM, mov, ttl)
% Docked viewer (dark): rows = classes; voltage movie (-dF), its ROI trace with the time cursor,
% glutamate movie (dF); gray colormap; slider = time
fig = figure('Name', ['STA movies ' ttl], 'NumberTitle', 'off', 'WindowStyle', 'docked');
[hI, hC, ax] = trigMovieLayout(fig, SM, mov);
sl = uicontrol(fig, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.1 0.005 0.8 0.025], ...
    'Min', SM.tV(1)*1000, 'Max', SM.tV(end)*1000, 'Value', 0, 'SliderStep', [1 10] / (numel(SM.tV) - 1));
addlistener(sl, 'Value', 'PostSet', @(~,~) show(sl.Value));
show(0);
    function show(tms)
        [~, kV] = min(abs(SM.tV*1000 - tms));
        trigMovieFrame(SM, hI, hC, ax, kV);
        drawnow limitrate
    end
end

function [hI, hC, ax] = trigMovieLayout(fig, SM, mov)
% 3 x 3 tiles: voltage movie | voltage ROI trace | glutamate movie, one row per class (dark, gray)
try, theme(fig, 'dark'); catch, set(fig, 'Color', 'k'); end
tl = tiledlayout(fig, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');   % left half voltage, right half glutamate
[clV, clG] = movieClims(SM, mov);
yl = [min(cellfun(@min, SM.voltTrace)) max(cellfun(@max, SM.voltTrace))]; yl = yl + [-0.05 0.05]*diff(yl);
ax = gobjects(3, 3); hI = gobjects(3, 2); hC = gobjects(3, 1);
for c = 1:3
    ax(c,1) = nexttile(tl, 4*(c-1) + 1); hI(c,1) = imagesc(ax(c,1), SM.volt{c}(:,:,1), clV); axis(ax(c,1), 'image', 'off'); colormap(ax(c,1), gray);
    hold(ax(c,1), 'on'); plot(ax(c,1), SM.roiV([1:end 1],1), SM.roiV([1:end 1],2), '-', 'Color', [1 0.55 0]); hold(ax(c,1), 'off');
    ax(c,2) = nexttile(tl, 4*(c-1) + 2); plot(ax(c,2), SM.tV*1000, SM.voltTrace{c}, '-', 'Color', [0.9 0.9 0.9], 'LineWidth', 1); hold(ax(c,2), 'on');
    hC(c) = xline(ax(c,2), SM.tV(1)*1000, '-', 'Color', [1 0.3 0.3], 'LineWidth', 1.5); hold(ax(c,2), 'off');
    xlim(ax(c,2), SM.tV([1 end])*1000); ylim(ax(c,2), yl); box(ax(c,2), 'off');
    ylabel(ax(c,2), 'voltage -\DeltaF (ROI)'); if c == 3, xlabel(ax(c,2), 'time from spike (ms)'); end
    title(ax(c,2), sprintf('%s, n = %d', SM.names{c}, SM.nEvents(c)), 'FontSize', 9);
    ax(c,3) = nexttile(tl, 4*(c-1) + 3, [1 2]); hI(c,2) = imagesc(ax(c,3), SM.glu{c}(:,:,1), clG); axis(ax(c,3), 'image', 'off'); colormap(ax(c,3), gray);
    hold(ax(c,3), 'on'); plot(ax(c,3), SM.roiG([1:end 1],1), SM.roiG([1:end 1],2), '-', 'Color', [1 0.55 0]); hold(ax(c,3), 'off');
end
cb = colorbar(ax(1,1)); cb.Label.String = 'voltage -\DeltaF'; cb = colorbar(ax(1,3)); cb.Label.String = 'glutamate \DeltaF';
title(tl, '\DeltaF relative to each pixel''s temporal median', 'FontSize', 9);
end

function trigMovieFrame(SM, hI, hC, ax, kV)
% Show voltage frame kV and the glutamate at its time in every row
for c = 1:3
    hI(c,1).CData = SM.volt{c}(:,:,kV); hI(c,2).CData = gluAt(SM, c, SM.tV(kV));
    hC(c).Value = SM.tV(kV)*1000;
    title(ax(c,1), sprintf('voltage, %+.0f ms', SM.tV(kV)*1000), 'FontSize', 9);
    title(ax(c,3), sprintf('glutamate, %+.0f ms', SM.tV(kV)*1000), 'FontSize', 9);
end
end

function G = gluAt(SM, c, t)
% Glutamate STA frame at time t (s), linearly interpolated between the grid frames
k = find(SM.tG <= t, 1, 'last');
if isempty(k), G = SM.glu{c}(:,:,1); return; end
if k >= numel(SM.tG), G = SM.glu{c}(:,:,end); return; end
a = (t - SM.tG(k)) / (SM.tG(k+1) - SM.tG(k));
G = (1 - a) * SM.glu{c}(:,:,k) + a * SM.glu{c}(:,:,k+1);
end

function [clV, clG] = movieClims(SM, mov)
% Colour limits of each channel: percentiles over all classes and frames
v = cellfun(@(m) m(1:7:end)', SM.volt, 'UniformOutput', false); v = [v{:}];
g = cellfun(@(m) m(1:7:end)', SM.glu, 'UniformOutput', false); g = [g{:}];
clV = prctile(v(isfinite(v)), mov.pctV); clG = prctile(g(isfinite(g)), mov.pctG);
if ~(diff(clV) > 0), clV = [-1 1]; end
if ~(diff(clG) > 0), clG = [-1 1]; end
end

function trigMovieMP4(SM, mov, fn)
% MP4 of the viewer layout (dark), every mov.step_ms of the voltage frames
fig = figure('Visible', 'off', 'Position', [50 50 1500 900]);
[hI, hC, ax] = trigMovieLayout(fig, SM, mov);
vw = VideoWriter(fn, 'MPEG-4'); vw.FrameRate = 20; open(vw);
step = max(1, round(mov.step_ms / 1000 / median(diff(SM.tV))));
for kV = 1:step:numel(SM.tV)
    trigMovieFrame(SM, hI, hC, ax, kV);
    drawnow; writeVideo(vw, getframe(fig));
end
close(vw); close(fig);
fprintf('  movie saved: %s\n', fn);
end

function p = gluLagOf(fp, sta)
% Section 3-0 lag line of recording fp (computed now, with a warning, when it is not in sta.gluLag)
k = [];
if isfield(sta, 'gluLag') && ~isempty(sta.gluLag), k = find(strcmp({sta.gluLag.fpath}, fp), 1); end
if ~isempty(k), p = sta.gluLag(k).p; return; end
if ~isfield(sta, 'lagP'), error('No glutamate lag for %s: run Section 3-0', fp); end
warning('gluLagOf:missing', '%s not in gluLag (Section 3-0): lag computed now', fp);
p = gluLagFit(fp, sta.lagP);
end

function [p, info] = gluLagFit(fp, P)
% Lag (ms) of the glutamate exposure start after Glu t_ax as a line in time, lag(t) = polyval(p, t):
% per-tile lags from the motion of the two cameras (motionLag), tiles with a correlation peak
% >= P.minR and within 4 ms of their median; fewer than 2 such tiles: the whole-recording lag
% (constant) if its peak >= P.minR, else no correction (warning). P: Section 3-0 parameters (lagP).
V = load(fullfile(fp, 'Volt_Result.mat'), 'Result'); V = V.Result;
Gm = loadRetry(fullfile(fp, 'Glu_Result.mat'), 'mc', 't_ax', 'tileLen');
L = motionLag(V.mc, V.t_ax, Gm.mc, Gm.t_ax, Gm.tileLen, P.A, struct('tau_ms', P.tau_ms, 'hp_s', P.hp_s));
e = cumsum(Gm.tileLen(:))'; s0 = [1, e(1:end-1) + 1]; e = min(e, numel(Gm.t_ax));
tMid = (Gm.t_ax(s0) + Gm.t_ax(e)) / 2; tMid = tMid(:)';
good = L.tileR >= P.minR & abs(L.tileLag_ms - median(L.tileLag_ms(L.tileR >= P.minR))) <= 4;
if nnz(good) >= 2
    p = polyfit(tMid(good), L.tileLag_ms(good), 1); how = sprintf('line over %d / %d tiles', nnz(good), numel(good));
elseif L.rPeak >= P.minR
    p = [0 L.lag_ms]; how = 'constant (whole recording)';
else
    p = [0 0]; how = 'none (no clear motion correlation)';
    warning('gluLagFit:none', '%s: no clear voltage-glutamate motion correlation (r %.2f); glutamate time not corrected', fp, L.rPeak);
end
info = struct('lag_ms', L.lag_ms, 'rPeak', L.rPeak, 'tileLag_ms', L.tileLag_ms, 'tileR', L.tileR, 'tileMid_s', tMid, ...
    'good', good, 'p', p, 'how', how);
end

function L = motionLag(mcV, tV, mcG, tG, tileLen, A, lag)
% Lag between the voltage and glutamate time axes from the rigid motion of the two cameras (Section 3-0)
rig = @(m) [mean(m(:, 1:end/2), 2), mean(m(:, end/2+1:end), 2)];   % [x patches, y patches] -> rigid [x y]
sV = rig(double(mcV)) * A';                              % voltage shift in the glutamate camera axes
sG = rig(double(mcG));
nV = min(size(sV,1), numel(tV)); tV = tV(1:nV); tV = tV(:); sV = sV(1:nV,:);
nG = min(size(sG,1), numel(tG)); tG = tG(1:nG); tG = tG(:); sG = sG(1:nG,:);
dtV = median(diff(tV)); dtG = median(diff(tG));
tE = [tV - dtV/2; tV(end) + dtV/2];                      % voltage frame edges
CI = [zeros(1,2); cumsum(sV * dtV, 1)];                  % integral of the voltage shift
hpG = @(x) x - movmean(x, max(3, round(lag.hp_s / dtG)), 1);
yG = hpG(sG); yG = (yG - mean(yG)) ./ std(yG);
r = nan(size(lag.tau_ms));
ok0 = []; P = cell(size(lag.tau_ms));
for k = 1:numel(lag.tau_ms)
    a = tG + lag.tau_ms(k)/1000; b = a + dtG;
    ok = a >= tE(1) & b <= tE(end);
    pv = nan(nG, 2);
    pv(ok,:) = (interp1(tE, CI, b(ok)) - interp1(tE, CI, a(ok))) / dtG;   % voltage shift over the exposure
    if isempty(ok0), ok0 = ok; else, ok0 = ok0 & ok; end
    P{k} = pv;
end
for k = 1:numel(lag.tau_ms)                              % same frames for every lag
    x = hpG(fillmissing(P{k}, 'nearest')); x = x(ok0,:); y = yG(ok0,:);
    x = (x - mean(x)) ./ std(x);
    r(k) = mean(sum(x .* y, 1) / (size(x,1) - 1));      % mean of the x and y correlations
end
[L.rPeak, ip] = max(r); L.lag_ms = lag.tau_ms(ip);
L.r0 = r(lag.tau_ms == 0); if isempty(L.r0), L.r0 = NaN; end
L.r = r;
% per glutamate tile
e = cumsum(tileLen(:))'; s0 = [1, e(1:end-1) + 1]; e = min(e, nG);
L.tileLag_ms = nan(1, numel(e)); L.tileR = nan(1, numel(e));
for t = 1:numel(e)
    in = false(nG, 1); in(s0(t):e(t)) = true; in = in & ok0;
    if nnz(in) < 100, continue; end
    rt = nan(size(r));
    for k = 1:numel(lag.tau_ms)
        x = hpG(fillmissing(P{k}, 'nearest')); x = x(in,:); y = yG(in,:);
        x = (x - mean(x)) ./ std(x); y = (y - mean(y)) ./ std(y);
        rt(k) = mean(sum(x .* y, 1) / (size(x,1) - 1));
    end
    [L.tileR(t), it] = max(rt); L.tileLag_ms(t) = lag.tau_ms(it);
end
end
