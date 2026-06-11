function VoltronResult = extractVoltronST(mov2_sub, mov2_mc_filt, opts, S_glu_in)
% extractVoltronST
%   Variant of extractGluSNFR3 designed for movies with blinking/large cells.
%
% KEY DIFFERENCES from extractGluSNFR3
%   1. Larger footprint radius: r_max default = 20 px (vs 8)
%   2. No NMF split — all blobs kept as single events
%   3. Detection uses temporal high-pass movie (mov - movmedian(mov,hpWin,3))
%      to find peaks and footprints that exceed threshold
%   4. Trace extraction (pseudo-inverse) uses original mov2_sub, not
%      the high-pass filtered movie
%
% PIPELINE
%   0. High-pass filter in time, then robust per-pixel z-score (MAD-based)
%   1. Seed pixels: 3x3-averaged z-score traces, sampled at seedStride
%   2. Spatial blob segmentation at each seed frame (no NMF split)
%   3. Centroid clustering  (r_min graph-connect + k-means for large clusters)
%   4. Build S_glu from mean accumulated z-score footprints, L2-normalised
%   5. Boundary filter + interactive synapse filter UI
%   6. Least-squares trace extraction  (S \ mov2_sub, original movie)
%   7. Diagnostic figures
%
% MODES
%   'detect'  (default) — full pipeline
%   'project' — project existing S_glu_in onto movie, skip detection
%
% USAGE
%   R = extractVoltronST(mov2_sub, mov2_mc_filt)
%   R = extractVoltronST(mov2_sub, mov2_mc_filt, opts)
%   R = extractVoltronST(mov2_sub, mov2_mc_filt, opts, S_glu_in)
%
% PARAMETERS (opts struct — all fields are optional, defaults shown)
%
%  --- Detection ---
%   .hpWin         movmedian window for temporal high-pass filter [50]
%                  Removes slow baseline/blinking. Increase for slower signals.
%   .smoothSigma   Gaussian spatial pre-filter sigma, px [1]
%                  Applied before z-scoring. Increase to reduce pixel noise.
%   .zSeedMin      Min z-score of the 3x3-averaged trace to be considered a
%                  seed peak [5.5]. Raise to reduce false positives.
%   .zThrFoot      Min z-score to include a pixel in a blob footprint [2.5]
%                  Lower grows footprint boundaries; raise to shrink them.
%   .minPeakDist_fr Min number of frames between seed peaks at the same pixel [10]
%   .minPeakProm   findpeaks MinPeakProminence on z-score trace [0.75]
%   .minArea       Min blob size to keep, pixels [2]
%   .seedStride    Spatial stride for the seed pixel scan [2]
%                  Stride 2 samples every other pixel (4x fewer findpeaks calls).
%                  Increase to 3–4 for large cells / fast runs.
%
%  --- Clustering & footprints ---
%   .r_max         Clustering radius: two events within r_max px of each
%                  other are merged into the same synapse. Also used as the
%                  radius around the synapse centroid within which blob pixels
%                  are accepted into the footprint, px [20]
%                  This is the main parameter for controlling footprint size.
%   .bound         Border margin excluded from synapse detection, px [10]
%   .minEventsSyn  Min number of events a synapse must have to be kept [1]
%   .footSigma     Gaussian blur applied to each footprint before L2-norm [0.75]
%                  Set to 0 to disable smoothing.
%   .usePeakWeight Weight each event's z-score by its peak value when
%                  accumulating the footprint [true]
%
%  --- Synapse filter (interactive or replay) ---
%   .synFilter     Saved filter struct with fields ampRange, countRange,
%                  areaRange. Applied non-interactively when doPlot=false. []
%   .saveSynFilter Save the interactively chosen filter to synFilterFile [false]
%   .synFilterFile Path for saving/loading synFilter struct ['']
%
%  --- Output & figures ---
%   .doPlot        Show all diagnostic figures and interactive UIs [true]
%   .uiSeedThresh  Interactive slider to adjust zSeedMin [true]
%   .uiClusterSize Interactive slider to adjust r_max [true]
%   .uiSynFilter   Interactive amplitude/count/area synapse filter [true]
%   .diagSaveDir   Directory for saved PNG diagnostics [pwd]
%   .diagPrefix    Filename prefix for saved PNGs ['VoltronDiag']
%   .diagZoom_pad  Pixels around centroid in zoomed footprint panels [15]
%   .diagMaxSyn    Max synapses shown in per-synapse trace figure [16]

%% ===================================================================
%  DEFAULTS
%% ===================================================================
if nargin < 3 || isempty(opts), opts = struct; end

    function v = go(f, d)
        if isfield(opts,f) && ~isempty(opts.(f)); v = opts.(f); else; v = d; end
    end

% Detection
opts.hpWin           = go('hpWin',           30);
opts.smoothSigma     = go('smoothSigma',     1);
opts.zSeedMin        = go('zSeedMin',        5.5);
opts.zThrFoot        = go('zThrFoot',        2.5);
opts.minPeakDist_fr  = go('minPeakDist_fr',  10);
opts.minPeakProm     = go('minPeakProm',     0.75);
opts.minArea         = go('minArea',         15);
opts.seedStride      = go('seedStride',      3);
% Clustering & footprints
opts.r_max           = go('r_max',           20);
opts.bound           = go('bound',           10);
opts.minEventsSyn    = go('minEventsSyn',    1);
opts.footSigma       = go('footSigma',       1);
opts.usePeakWeight   = go('usePeakWeight',   true);
% Synapse filter
opts.synFilter       = go('synFilter',       []);
opts.saveSynFilter   = go('saveSynFilter',   false);
opts.synFilterFile   = go('synFilterFile',   '');
% Output & figures
opts.doPlot          = go('doPlot',          true);
opts.uiSeedThresh    = go('uiSeedThresh',    true);
opts.uiClusterSize   = go('uiClusterSize',   true);
opts.uiSynFilter     = go('uiSynFilter',     true);
opts.diagSaveDir     = go('diagSaveDir',     pwd);
opts.diagPrefix      = go('diagPrefix',      'VoltronDiag');
opts.diagZoom_pad    = go('diagZoom_pad',    15);
opts.diagMaxSyn      = go('diagMaxSyn',      16);

if ~exist(opts.diagSaveDir,'dir'), mkdir(opts.diagSaveDir); end

%% ===================================================================
%  PROJECT MODE
%% ===================================================================
if isfield(opts,'mode') && strcmpi(opts.mode,'project')
    if nargin < 4 || isempty(S_glu_in)
        error('extractVoltronST: project mode requires S_glu_in.');
    end
    mov = double(mov2_sub);
    [H,W,T] = size(mov);
    [Hs,Ws,Nsyn] = size(S_glu_in);
    assert(H==Hs && W==Ws,'S_glu_in spatial dims must match movie.');
    Smat = reshape(double(S_glu_in), H*W, Nsyn);
    Mmat = reshape(mov, H*W, T);
    Mmat(~isfinite(Mmat)) = 0;
    T_glu = Smat \ Mmat;
    F0_glu = []; dFF_glu = [];
    if ~isempty(mov2_mc_filt)
        M0 = reshape(double(mov2_mc_filt), H*W, T);
        M0(~isfinite(M0)) = 0;
        F0_glu  = Smat \ M0;
        dFF_glu = T_glu ./ mean(F0_glu, 2);
    end
    VoltronResult = struct('S_glu',S_glu_in,'T_glu',T_glu,...
        'F0_glu',F0_glu,'dFF_glu',dFF_glu,'opts',opts);
    return
end

%% ===================================================================
%  STEP 0: PREPROCESS — high-pass filter then robust z-score
%% ===================================================================
mov = double(mov2_sub);
[H,W,T] = size(mov);
minPeakDist = max(1, round(opts.minPeakDist_fr));

% High-pass filter in time: remove slow baseline / blinking artefacts
fprintf('Step 0: computing high-pass filtered movie (movmedian window=%d)...\n', opts.hpWin);
mov_hp = mov - movmedian(mov, opts.hpWin, 3);

mov_hp_sm = imgaussfilt(mov_hp, opts.smoothSigma);

% Robust per-pixel z-score on high-pass movie
medImg   = median(mov_hp_sm, 3);
madImg   = mad(mov_hp_sm, 1, 3);
sigmaImg = 1.4826 * madImg + eps;
zMov     = (mov_hp_sm - medImg) ./ sigmaImg;

% Max-projection for UI
zMaxProj = max(zMov, [], 3, 'omitnan');
nzv = zMaxProj(zMaxProj ~= 0);
if ~isempty(nzv), zMaxProj(zMaxProj==0) = median(nzv); end

% 3x3 spatial average for peak seeding: smoother traces, fewer false peaks
zMov3 = convn(zMov, ones(3,3,1)/9, 'same');
% Pre-compute temporal max so the seed loop can skip pixels cheaply
zMax3 = max(zMov3, [], 3);

%% ===================================================================
%  UI: SEED THRESHOLD
%% ===================================================================
if opts.doPlot && opts.uiSeedThresh
    zMin = 0;
    zMax = max(6, ceil(max(zMov(:),[],'omitnan')));

    fig = uifigure('Name','Adjust zSeedMin','Position',[80 80 920 590]);
    ax1 = uiaxes(fig,'Position',[20  140 430 415]);
    ax2 = uiaxes(fig,'Position',[465 140 430 415]);
    sld = uislider(fig,'Position',[20 95 880 3],...
                   'Limits',[zMin zMax],'Value',opts.zSeedMin);
    lbl = uilabel(fig,'Position',[20 62 600 26],...
                  'Text',sprintf('zSeedMin = %.2f', sld.Value),'FontSize',13);
    btn = uibutton(fig,'Text','Accept','Position',[720 18 180 36],...
                   'FontSize',13,'ButtonPushedFcn',@(~,~) uiresume(fig)); %#ok<NASGU>
    sld.ValueChangedFcn = @(src,~) redrawSeed(src);
    redrawSeed(sld);
    uiwait(fig);
    opts.zSeedMin = sld.Value;
    close(fig);
    fprintf('zSeedMin = %.2f\n', opts.zSeedMin);
end

%% ===================================================================
%  STEP 1: PIXEL-PEAK SEEDS  (on high-pass z-score movie)
%% ===================================================================
stride = opts.seedStride;
fprintf('Step 1: seeding peaks (3x3 avg z-score, stride=%d)...\n', stride);
seedByFrame = cell(T,1);
for y = 1:stride:H
    for x = 1:stride:W
        if zMax3(y,x) < opts.zSeedMin, continue; end   % cheap pre-check
        tr = squeeze(zMov3(y,x,:));
        [pks,locs] = findpeaks(tr,...
            'MinPeakDistance',   minPeakDist,...
            'MinPeakProminence', opts.minPeakProm);
        if isempty(locs), continue; end
        locs = locs(pks >= opts.zSeedMin);
        if isempty(locs), continue; end
        lin = sub2ind([H,W],y,x);
        for k = 1:numel(locs)
            seedByFrame{locs(k)}(end+1) = lin; %#ok<AGROW>
        end
    end
end
fprintf('  %d seed pixels found.\n', sum(cellfun(@numel,seedByFrame)));

%% ===================================================================
%  STEP 2: BLOB SEGMENTATION  (no NMF split)
%% ===================================================================
% Each detected blob is kept as a single event regardless of size.
% The high-pass z-score movie is used for thresholding.
fprintf('Step 2: segmenting blobs...\n');

eventList = struct('frame',{},'pixIdx',{},'centroid_xy',{},...
                   'peakZ',{},'pixVal',{},'meanZ',{});

compCount = 0;

for t = 1:T
    seeds = seedByFrame{t};
    if isempty(seeds), continue; end

    zFrame = zMov(:,:,t);
    bw = zFrame > opts.zThrFoot;
    bw = bwareaopen(bw, opts.minArea);
    bw = imopen(bw, strel('disk',1));
    bw = imfill(bw,'holes');
    CC = bwconncomp(bw,8);
    if CC.NumObjects == 0, continue; end

    L = labelmatrix(CC);
    seedLabels = unique(L(seeds));
    seedLabels(seedLabels==0) = [];
    if isempty(seedLabels), continue; end

    for lab = seedLabels(:)'
        pixIdx = find(L==lab);
        if numel(pixIdx) < opts.minArea, continue; end
        [r,c]  = ind2sub([H,W], pixIdx);
        pixVal = zFrame(pixIdx);
        compCount = compCount + 1;
        eventList(compCount) = buildEvent(t, pixIdx, c, r, pixVal);
    end
end

fprintf('  %d events extracted.\n', numel(eventList));

if isempty(eventList)
    warning('No events detected. Returning empty result.');
    VoltronResult = struct('eventList',[],'S_glu',[],'T_glu',[],...
                           'dFF_glu',[],'opts',opts);
    return
end

%% ===================================================================
%  STEP 3: CLUSTER EVENT CENTROIDS → SYNAPSES
%% ===================================================================
fprintf('Step 3: clustering %d event centroids...\n', numel(eventList));
cent   = reshape([eventList.centroid_xy],2,[])';   % Nevent x 2
Nevent = size(cent,1);

% Optional UI to tune r_max
if opts.doPlot && opts.uiClusterSize && Nevent > 0
    % Background image: mean of mov2_mc_filt, or mean of mov2_sub as fallback
    if ~isempty(mov2_mc_filt)
        bgImgCluster = mean(double(mov2_mc_filt), 3);
    else
        bgImgCluster = mean(double(mov2_sub), 3);
    end

    fig = uifigure('Name','Adjust r\_max — clustering radius (events within r\_max are merged)',...
                   'Position',[80 80 920 590]);
    ax  = uiaxes(fig,'Position',[20 115 880 445]);
    axis(ax,'ij'); axis(ax,'equal'); xlim(ax,[1 W]); ylim(ax,[1 H]);
    imagesc(ax, bgImgCluster); colormap(ax, gray(256));
    hold(ax,'on');
    b = opts.bound;
    rectangle(ax,'Position',[b+1,b+1,W-2*b,H-2*b],...
        'EdgeColor','y','LineWidth',1.5,'LineStyle','--');
    sld = uislider(fig,'Position',[20 78 880 3],...
                   'Limits',[0.5 max(30, opts.r_max*5)],'Value',opts.r_max);
    lbl = uilabel(fig,'Position',[20 50 600 26],...
                  'Text',sprintf('r_max = %.2f px',opts.r_max),'FontSize',13);
    btn = uibutton(fig,'Text','Accept','Position',[720 14 180 36],...
                   'FontSize',13,'ButtonPushedFcn',@(~,~) uiresume(fig)); %#ok<NASGU>
    fig.UserData.ax = ax;
    fig.UserData.sh = [];
    sld.ValueChangedFcn = @(src,~) redrawCluster(fig,lbl,src.Value);
    redrawCluster(fig,lbl,sld.Value);
    uiwait(fig);
    opts.r_max = sld.Value;
    close(fig);
    fprintf('r_max = %.2f px\n', opts.r_max);
end

% Pairwise distance clustering: events within r_max px are linked (single-linkage)
A        = pdist2(cent, cent, 'squaredeuclidean') <= opts.r_max^2;
clusters = accumarray(conncomp(graph(A))',(1:Nevent)',[],@(v){v});

% Split any cluster whose radius (max dist from centroid) exceeds 1.5 * r_max
synID      = zeros(Nevent,1);
synCounter = 0;
for c = 1:numel(clusters)
    idx = clusters{c};
    pts = cent(idx,:);
    ctr = mean(pts, 1);
    radius = sqrt(max(sum((pts - ctr).^2, 2)));

    if radius <= 1.2 * opts.r_max || numel(idx) <= 2
        synCounter = synCounter + 1;
        synID(idx) = synCounter;
    else
        K = max(2, round(radius / opts.r_max));
        try
            lab = kmeans(pts, K, 'Replicates', 5, 'MaxIter', 300);
        catch
            synCounter = synCounter + 1;
            synID(idx) = synCounter;
            continue
        end
        for k = 1:K
            si = idx(lab == k);
            if isempty(si), continue; end
            synCounter = synCounter + 1;
            synID(si) = synCounter;
        end
    end
end

fprintf('  %d events → %d synapses (r_max=%.1f px)\n',...
    Nevent, synCounter, opts.r_max);

if opts.doPlot && synCounter > 0
    fh = figure('Name','Event centroids by synapse');
    hold on; axis ij; axis equal; xlim([1 W]); ylim([1 H]);
    b = opts.bound;
    rectangle('Position',[b+1,b+1,W-2*b,H-2*b],...
        'EdgeColor','w','LineWidth',2,'LineStyle','--');
    cmap_s = hsv(synCounter);
    for s = 1:synCounter
        ei = find(synID==s); if isempty(ei), continue; end
        scatter(cent(ei,1),cent(ei,2),14,cmap_s(s,:),'filled');
    end
    title('Event centroids colored by synapse ID');
    set(gca,'Color','k'); hold off;
    saveDiag(fh,opts,'centroids_by_synapse');
end

%% ===================================================================
%  STEP 4: BUILD S_glu — mean accumulated z-score footprints
%% ===================================================================
% Footprint pixels are cropped to r_max (=15 px) of each synapse centroid.
% z-scores come from the high-pass filtered movie.
fprintf('Step 4: building synapse footprints (r_max crop = %.1f px)...\n', opts.r_max);

% Pre-compute synapse centroids
synCentroid0 = nan(synCounter, 2);
for s = 1:synCounter
    ei = find(synID == s);
    if ~isempty(ei)
        synCentroid0(s,:) = mean(cent(ei,:), 1);
    end
end

S_glu0       = zeros(H,W,synCounter,'double');
synEventCnt0 = zeros(synCounter,1);

for e = 1:Nevent
    s = synID(e);
    if s == 0, continue; end
    synEventCnt0(s) = synEventCnt0(s) + 1;

    pixIdx_d = double(eventList(e).pixIdx(:));
    pixVal_d = double(eventList(e).pixVal(:));

    % Restrict footprint to 1.2 * r_max from synapse centroid
    cx = synCentroid0(s,1);
    cy = synCentroid0(s,2);
    [pr, pc] = ind2sub([H,W], pixIdx_d);
    inRadius = (pc - cx).^2 + (pr - cy).^2 <= (1.2 * opts.r_max)^2;
    pixIdx_d = pixIdx_d(inRadius);
    pixVal_d = pixVal_d(inRadius);
    if isempty(pixIdx_d), continue; end

    if opts.usePeakWeight
        w = max(eventList(e).peakZ, 0);
        if ~isfinite(w) || w==0, w=1; end
    else
        w = 1;
    end

    fp_slice = zeros(H,W);
    fp_slice(pixIdx_d) = w * pixVal_d;
    S_glu0(:,:,s) = S_glu0(:,:,s) + fp_slice;
end

% Divide by event count → mean weighted z-score footprint
for s = 1:synCounter
    if synEventCnt0(s) > 0
        S_glu0(:,:,s) = S_glu0(:,:,s) / synEventCnt0(s);
    end
end

% Optional spatial smoothing
if opts.footSigma > 0
    for s = 1:synCounter
        S_glu0(:,:,s) = imgaussfilt(S_glu0(:,:,s), opts.footSigma);
    end
end

% L2 normalise
for s = 1:synCounter
    v = S_glu0(:,:,s);
    S_glu0(:,:,s) = v / (sqrt(sum(v(:).^2)) + eps);
end

%% ===================================================================
%  STEP 5: BOUNDARY FILTER + INTERACTIVE SYNAPSE FILTER
%% ===================================================================
fprintf('Step 5: boundary filter + synapse selection...\n');
bd   = opts.bound;
cSum = squeeze(sum(S_glu0(bd+1:H-bd, bd+1:W-bd, :), [1 2]));
keep = cSum > 0;

S_glu   = S_glu0(:,:,keep);
newNsyn = size(S_glu,3);

mapOld = zeros(synCounter,1);
mapOld(find(keep)) = 1:newNsyn;
newSynID = mapOld(synID);

% Recompute centroid / event count for kept synapses
synEventCount = zeros(newNsyn,1);
synCentroid   = nan(newNsyn,2);
for s = 1:newNsyn
    ei = find(newSynID==s);
    synEventCount(s) = numel(ei);
    if ~isempty(ei), synCentroid(s,:) = mean(cent(ei,:),1); end
end

% Footprint features (vectorized over synapses)
Sflat   = reshape(S_glu, H*W, newNsyn);   % pixels x synapses
posMask = Sflat > 0;
synArea = sum(posMask, 1)';
synAmp  = sum(Sflat .* posMask, 1)' ./ max(synArea, 1);

% Amplitude proxy: dot high-pass zMov onto each footprint
amp  = max(reshape(zMov,H*W,T)' * reshape(S_glu,H*W,newNsyn), [], 1)';
cnt  = synEventCount;
area = synArea;

% Apply saved filter without UI
keepUI = true(newNsyn,1);
if ~isempty(opts.synFilter) && ~(opts.doPlot && opts.uiSynFilter)
    fr = opts.synFilter;
    keepUI = (amp  >= fr.ampRange(1)   & amp  <= fr.ampRange(2))  & ...
             (cnt  >= fr.countRange(1) & cnt  <= fr.countRange(2)) & ...
             (area >= fr.areaRange(1)  & area <= fr.areaRange(2));
end

% Background image for the filter UI: temporal mean of mov2_mc_filt.
% Falls back to mean(mov2_sub) if mov2_mc_filt was not supplied.
if ~isempty(mov2_mc_filt)
    bgImgFilter = mean(double(mov2_mc_filt), 3);
else
    bgImgFilter = mean(double(mov2_sub), 3);
end

% Interactive filter UI
if opts.doPlot && opts.uiSynFilter && newNsyn > 0
    fig = uifigure('Name','Synapse filter','Position',[80 80 940 590]);
    ax  = uiaxes(fig,'Position',[20 235 900 330]);
    axis(ax,'ij'); axis(ax,'equal'); xlim(ax,[1 W]); ylim(ax,[1 H]);
    title(ax,'Red = kept     Green = filtered  |  background: mean(mov2\_mc\_filt)');
    imagesc(ax, bgImgFilter); colormap(ax,gray(256)); hold(ax,'on');
    b = opts.bound;
    rectangle(ax,'Position',[b+1,b+1,W-2*b,H-2*b],...
        'EdgeColor','y','LineWidth',1.5,'LineStyle','--');
    hKeep = scatter(ax,synCentroid(:,1),synCentroid(:,2),28,'o',...
        'MarkerEdgeColor',[1 0 0],'MarkerFaceColor','none','LineWidth',1.5);
    hDrop = scatter(ax,synCentroid(:,1),synCentroid(:,2),28,'o',...
        'MarkerEdgeColor',[0 1 0],'MarkerFaceColor','none','LineWidth',1.5);
    hold(ax,'off');

    ampL  = prcLims(amp,  0,100,eps);
    cntL  = prcLims(cnt,  0,100,1);
    areaL = prcLims(area, 0,100,1);

    [sAmpLo, sAmpHi]  = sliderPair(fig,195,'Amplitude',  ampL,  ampL);
    [sCntLo, sCntHi]  = sliderPair(fig,148,'Event count', cntL,  cntL);
    [sAreaLo,sAreaHi] = sliderPair(fig,100,'Area (px)',   areaL, areaL);
    lbl = uilabel(fig,'Position',[20 60 680 26],'Text','','FontSize',12);
    btn = uibutton(fig,'Text','Accept','Position',[730 15 190 38],...
                   'FontSize',13,'ButtonPushedFcn',@(~,~) uiresume(fig)); %#ok<NASGU>

    for sh = [sAmpLo sAmpHi sCntLo sCntHi sAreaLo sAreaHi]
        sh.ValueChangedFcn = @(~,~) redrawFilter();
    end
    redrawFilter();
    uiwait(fig);

    a1=sAmpLo.Value; a2=sAmpHi.Value;
    c1=sCntLo.Value; c2=sCntHi.Value;
    r1=sAreaLo.Value; r2=sAreaHi.Value;
    opts.synFilter = struct('ampRange',sort([a1 a2]),...
                            'countRange',sort([c1 c2]),...
                            'areaRange',sort([r1 r2]));
    if opts.saveSynFilter && ~isempty(opts.synFilterFile)
        synFilter = opts.synFilter; %#ok<NASGU>
        save(opts.synFilterFile,'synFilter');
    end
    keepUI = (amp>=min(a1,a2) & amp<=max(a1,a2)) & ...
             (cnt>=min(c1,c2) & cnt<=max(c1,c2)) & ...
             (area>=min(r1,r2) & area<=max(r1,r2));
    close(fig);
end

% Apply UI mask
S_glu         = S_glu(:,:,keepUI);
synCentroid   = synCentroid(keepUI,:);
synEventCount = synEventCount(keepUI);
newNsyn       = size(S_glu,3);
mapUI         = zeros(numel(keepUI),1);
mapUI(find(keepUI)) = 1:newNsyn;
valid = newSynID>0;
tmp   = zeros(size(newSynID));
tmp(valid) = mapUI(newSynID(valid));
newSynID = tmp;

% Hard minimum-events cut
keepCnt       = synEventCount >= opts.minEventsSyn;
fprintf('  Removing %d synapses with < %d events.\n', sum(~keepCnt), opts.minEventsSyn);
S_glu         = S_glu(:,:,keepCnt);
synCentroid   = synCentroid(keepCnt,:);
synEventCount = synEventCount(keepCnt);
newNsyn       = size(S_glu,3);
map2  = zeros(numel(keepCnt),1);
map2(find(keepCnt)) = 1:newNsyn;
tmp   = zeros(size(newSynID));
valid = newSynID>0;
tmp(valid) = map2(newSynID(valid));
newSynID = tmp;

synEventRate = synEventCount / max(T, 1);   % events per frame
fprintf('  Kept %d synapses.  Median rate: %.4f ev/fr\n', newNsyn, median(synEventRate));

%% ===================================================================
%  STEP 6: OVERVIEW VISUALISATION
%% ===================================================================
Smap_counts = zeros(H,W,'single');
for e = 1:Nevent
    sn = newSynID(e); if sn==0, continue; end
    Smap_counts(eventList(e).pixIdx) = Smap_counts(eventList(e).pixIdx) + 1;
end

RGB = [];
if opts.doPlot && newNsyn > 0
    fh = figure('Name','Event count map');
    imagesc(Smap_counts); axis image; axis off; colormap hot; colorbar;
    title('Event count map — kept synapses');
    saveDiag(fh,opts,'event_count_map');

    fh = figure('Name','Synapse footprints');
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    show_footprnt(S_glu, mean(mov2_sub,3));
    sgtitle('Synapse footprints (turbo) | overlay on mean movie');
    saveDiag(fh,opts,'synapse_footprints');

    fh = figure('Name','Footprints colored by event rate');
    tiledlayout(1,1);
    show_footprint_heatmap(S_glu, synEventRate);
    title('Footprints colored by event rate (ev/frame)');
    saveDiag(fh,opts,'synapse_footprints_eventrate');

    cmap_rgb = hsv(newNsyn);
    RGB = zeros(H,W,3,'double');
    for s = 1:newNsyn
        fp = S_glu(:,:,s);
        if max(fp(:))==0, continue; end
        V  = fp/(max(fp(:))+eps);
        RGB = max(RGB, bsxfun(@times,reshape(cmap_rgb(s,:),1,1,3),repmat(V,[1 1 3])));
    end
end

%% ===================================================================
%  STEP 7: COMPUTE T_glu  (pseudo-inverse on original mov2_sub)
%% ===================================================================
% NOTE: footprints (S_glu) were found using the high-pass filtered movie,
% but trace extraction uses the original mov2_sub to preserve slow signal.
fprintf('Step 7: computing T_glu via pseudo-inverse on original mov2_sub...\n');
usedPix = find(tovec(max(S_glu,[],3)) > 0);
Smat    = tovec(S_glu);     Smat = Smat(usedPix,:);   % Npix x Nsyn
Msub    = tovec(double(mov2_sub)); Msub = Msub(usedPix,:);
Msub(~isfinite(Msub)) = 0;
T_glu = Smat \ Msub;                                  % Nsyn x T

F0_glu = []; dFF_glu = [];
if ~isempty(mov2_mc_filt)
    M0 = tovec(double(mov2_mc_filt)); M0 = M0(usedPix,:);
    M0(~isfinite(M0)) = 0;
    F0_glu  = Smat \ M0;
    dFF_glu = T_glu ./ mean(F0_glu, 2);
end

%% ===================================================================
%  STEP 8: TRACE DIAGNOSTIC FIGURES
%% ===================================================================
if opts.doPlot && newNsyn > 0
    tvec  = 0:T-1;   % frame indices
    pad   = opts.diagZoom_pad;
    nShow = min(newNsyn, opts.diagMaxSyn);
    nCols = 4; nRows = ceil(nShow/2);

    fh = figure('Name','Per-synapse footprint + T\_glu trace',...
                'Position',[50 50 1400 nRows*210]);
    for s = 1:nShow
        xc_ = min(max(round(synCentroid(s,1)),1),W);
        yc_ = min(max(round(synCentroid(s,2)),1),H);
        xrng = max(1,xc_-pad):min(W,xc_+pad);
        yrng = max(1,yc_-pad):min(H,yc_+pad);
        row = ceil(s/2); col = mod(s-1,2)*2+1;

        ax1 = subplot(nRows,nCols,(row-1)*nCols+col);
        imagesc(ax1,S_glu(yrng,xrng,s)); axis(ax1,'image'); axis(ax1,'off');
        colormap(ax1,hot);
        hold(ax1,'on');
        plot(ax1,xc_-xrng(1)+1,yc_-yrng(1)+1,...
            'c+','MarkerSize',10,'LineWidth',1.5);
        hold(ax1,'off');
        title(ax1,sprintf('Syn %d',s),'FontSize',8);

        ax2 = subplot(nRows,nCols,(row-1)*nCols+col+1);
        plot(ax2,tvec,T_glu(s,:),'b','LineWidth',0.8);
        xlabel(ax2,'Frame','FontSize',7);
        title(ax2,sprintf('n=%d  %.4f ev/fr',...
            synEventCount(s),synEventRate(s)),'FontSize',8);
        set(ax2,'FontSize',7);
    end
    sgtitle('Per-synapse footprint (left) | T\_glu trace (right)','FontSize',11);
    saveDiag(fh,opts,'per_synapse_footprint_trace');

    nc2=ceil(sqrt(newNsyn)); nr2=ceil(newNsyn/nc2);
    fh=figure('Name','All synapse footprints',...
              'Position',[50 50 nc2*110 nr2*110]);
    for s = 1:newNsyn
        xc_=min(max(round(synCentroid(s,1)),1),W);
        yc_=min(max(round(synCentroid(s,2)),1),H);
        xrng=max(1,xc_-pad):min(W,xc_+pad);
        yrng=max(1,yc_-pad):min(H,yc_+pad);
        subplot(nr2,nc2,s);
        imagesc(S_glu(yrng,xrng,s)); axis image; axis off; colormap hot;
        title(sprintf('#%d n=%d',s,synEventCount(s)),'FontSize',6);
    end
    sgtitle('All footprints (zoomed, hot)','FontSize',10);
    saveDiag(fh,opts,'all_footprints_grid');

    fprintf('All figures saved to: %s\n', opts.diagSaveDir);
end

%% ===================================================================
%  OUTPUT STRUCT
%% ===================================================================
VoltronResult = struct();
VoltronResult.eventList       = eventList;
VoltronResult.Smap_event      = Smap_counts;
VoltronResult.S_glu           = S_glu;
VoltronResult.synCentroid     = synCentroid;
VoltronResult.synEventCount   = synEventCount;
VoltronResult.synEventRate    = synEventRate;
VoltronResult.T_glu           = T_glu;
VoltronResult.F0_glu          = F0_glu;
VoltronResult.dFF_glu         = dFF_glu;
VoltronResult.RGB             = RGB;
VoltronResult.opts            = opts;
VoltronResult.nFrames         = T;

%% ===================================================================
%  LOCAL HELPERS
%% ===================================================================

    function ev = buildEvent(t, pixIdx, c, r, pixVal)
        ev.frame       = t;
        ev.pixIdx      = pixIdx;
        ev.centroid_xy = [mean(c) mean(r)];
        ev.peakZ       = max(pixVal);
        ev.pixVal      = pixVal;
        ev.meanZ       = mean(pixVal);
    end

    function saveDiag(fh, o, tag)
        fname = fullfile(o.diagSaveDir,...
            sprintf('%s_%s.png', o.diagPrefix, tag));
        try
            exportgraphics(fh, fname, 'Resolution', 150);
            fprintf('  Saved: %s\n', fname);
        catch ME
            warning('saveDiag: %s', ME.message);
        end
    end

    function L = prcLims(v, pLo, pHi, epsPad)
        v = v(isfinite(v(:)));
        if isempty(v), L=[0 1]; return; end
        L = [prctile(v,pLo) prctile(v,pHi)];
        if L(1)==L(2), L=[L(1) L(1)+epsPad]; end
    end

    function redrawSeed(src)
        zthr = src.Value;
        lbl.Text = sprintf('zSeedMin = %.2f', zthr);
        imagesc(ax1, max(zMov>zthr,[],3)); colormap(ax1,gray);
        axis(ax1,'image'); axis(ax1,'off');
        title(ax1,sprintf('max(zMov > %.2f, [], 3)',zthr),'FontSize',9);
        imagesc(ax2, zMaxProj); colormap(ax2,turbo);
        axis(ax2,'image'); axis(ax2,'off');
        title(ax2,'Max z-score projection (high-pass)','FontSize',9);
        drawnow;
    end

    function redrawCluster(figH, lblH, rmax)
        lblH.Text = sprintf('r_max = %.2f px', rmax);
        centL    = cent;
        N_       = size(centL,1);
        axL      = figH.UserData.ax;

        % Pairwise distance clustering + radius split (mirrors Step 3)
        AL   = pdist2(centL, centL, 'squaredeuclidean') <= rmax^2;
        cls_ = accumarray(conncomp(graph(AL))',(1:N_)',[],@(v){v});
        synID_ = zeros(N_,1);
        sc_    = 0;
        for ci_ = 1:numel(cls_)
            ii_  = cls_{ci_};
            pts_ = centL(ii_,:);
            ctr_ = mean(pts_,1);
            rad_ = sqrt(max(sum((pts_-ctr_).^2,2)));
            if rad_ <= 1.5*rmax || numel(ii_) <= 2
                sc_ = sc_+1; synID_(ii_) = sc_;
            else
                K_ = max(2, round(rad_/rmax));
                try
                    lb_ = kmeans(pts_,K_,'Replicates',3,'MaxIter',200);
                catch
                    sc_ = sc_+1; synID_(ii_) = sc_; continue
                end
                for k_=1:K_
                    si_=ii_(lb_==k_); if isempty(si_), continue; end
                    sc_=sc_+1; synID_(si_)=sc_;
                end
            end
        end

        cm  = hsv(max(sc_,1));
        rng(0); perm = randperm(max(sc_,1));
        cm  = cm(perm,:);
        C   = cm(synID_,:);

        % Always redraw: background image + scatter on top
        cla(axL);
        imagesc(axL, bgImgCluster); colormap(axL, gray(256));
        hold(axL,'on');
        b_ = opts.bound;
        rectangle(axL,'Position',[b_+1,b_+1,W-2*b_,H-2*b_],...
            'EdgeColor','y','LineWidth',1.5,'LineStyle','--');
        h1 = scatter(axL, centL(:,1), centL(:,2), 14, C, 'filled', 'o');
        hold(axL,'off');
        figH.UserData.sh = h1;
        title(axL, sprintf('r_{max}=%.1f px | %d synapses', rmax, sc_), ...
              'FontSize',10);
        drawnow;
    end

    function redrawFilter()
        a1_=sAmpLo.Value;  a2_=sAmpHi.Value;
        c1_=sCntLo.Value;  c2_=sCntHi.Value;
        r1_=sAreaLo.Value; r2_=sAreaHi.Value;
        ku = (amp>=min(a1_,a2_) & amp<=max(a1_,a2_)) & ...
             (cnt>=min(c1_,c2_) & cnt<=max(c1_,c2_)) & ...
             (area>=min(r1_,r2_) & area<=max(r1_,r2_));
        hKeep.XData=synCentroid(ku,1);  hKeep.YData=synCentroid(ku,2);
        hDrop.XData=synCentroid(~ku,1); hDrop.YData=synCentroid(~ku,2);
        lbl.Text=sprintf('Kept: %d / %d synapses',sum(ku),numel(ku));
        drawnow limitrate;
    end

    function [sLo,sHi] = sliderPair(figH, yPos, label, lims, v0)
        uilabel(figH,'Text',label,...
            'Position',[15 yPos 130 22],'FontSize',11);
        uilabel(figH,'Text','min',...
            'Position',[150 yPos-18 35 20],'FontSize',9);
        sLo = uislider(figH,'Position',[150 yPos+6 320 3],...
                       'Limits',lims,'Value',v0(1));
        uilabel(figH,'Text','max',...
            'Position',[490 yPos-18 35 20],'FontSize',9);
        sHi = uislider(figH,'Position',[490 yPos+6 320 3],...
                       'Limits',lims,'Value',v0(2));
    end

end % extractVoltronST
