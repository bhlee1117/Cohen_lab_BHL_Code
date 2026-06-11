function GluResult = extractGluSNFR3(mov2_sub, mov2_mc_filt, exposuretime2, opts, S_glu_in)
% extractGluSNFR3
%   Detect iGluSNFR synaptic events, cluster into synapses, build S_glu,
%   and compute temporal traces T_glu by least-squares projection.
%
% PIPELINE
%   0. Robust per-pixel z-score movie (MAD-based)
%   1. Pixel-peak seeding  (findpeaks on z-score traces)
%   2. Spatial blob segmentation of each seed frame
%      → blobs larger than opts.eventMaxArea are NMF-split (pixels x time)
%        into K sub-events, each with its own pixIdx + centroid
%   3. Centroid clustering  (r_min graph-connect + k-means for large clusters)
%   4. Build S_glu from mean accumulated z-score footprints, L2-normalised
%   5. Boundary filter + interactive synapse filter UI
%   6. Least-squares trace extraction  (S \ movie)
%   7. Diagnostic figures (uses show_footprnt / show_footprint_heatmap)
%
% MODES
%   'detect'  (default) — full pipeline
%   'project' — project existing S_glu_in onto movie, skip detection
%
% USAGE
%   R = extractGluSNFR3(mov2_sub, mov2_mc_filt, exposuretime2)
%   R = extractGluSNFR3(mov2_sub, mov2_mc_filt, exposuretime2, opts)
%   R = extractGluSNFR3(mov2_sub, mov2_mc_filt, exposuretime2, opts, S_glu_in)
%
% KEY opts  (all have defaults — pass empty struct to use all defaults)
%   .smoothSigma     Gaussian pre-filter (px)               [1]
%   .zSeedMin        z-score peak threshold                  [3]
%   .zThrFoot        z-score footprint threshold             [2.5]
%   .minPeakDist_s   min peak separation (s)                [0.2]
%   .minPeakProm     findpeaks MinPeakProminence             [0.75]
%   .minArea         min blob area (px)                      [2]
%   .eventMaxArea    blob area that triggers NMF split (px)  [40]
%   .nmfReps         nnmf replicates                         [5]
%   .nmfMaxIter      nnmf max iterations                     [300]
%   .nmfMinPix       min px per NMF component to keep        [4]
%   .r_min           centroid graph radius (px)              [3]
%   .r_max           max cluster bbox before k-means (px)   [8]
%   .bound           boundary exclusion margin (px)          [10]
%   .minEventsSyn    min events per kept synapse             [3]
%   .footSigma       footprint Gaussian blur sigma (0=off)   [0.75]
%   .usePeakWeight   weight events by peakZ                  [true]
%   .doPlot          show all figures                        [true]
%   .uiSeedThresh    interactive zSeedMin slider             [true]
%   .uiClusterSize   interactive r_max slider                [true]
%   .uiSynFilter     interactive synapse filter              [true]
%   .diagSaveDir     folder for saved PNG diagnostics        [pwd]
%   .diagPrefix      filename prefix for saved PNGs          ['GluDiag']
%   .diagZoom_pad    px around centroid for zoomed panels    [15]
%   .diagMaxSyn      max synapses in per-synapse figure      [16]
%   .showNMFDiag     show NMF blob-split diagnostics         [true]
%   .nmfDiagMaxEvent max blob splits shown in diagnostic     [6]

%% ===================================================================
%  DEFAULTS
%% ===================================================================
if nargin < 4 || isempty(opts), opts = struct; end

    function v = go(f, d)
        if isfield(opts,f) && ~isempty(opts.(f)); v = opts.(f); else; v = d; end
    end

frameRate = 1 / exposuretime2;

opts.smoothSigma     = go('smoothSigma',     1);
opts.zSeedMin        = go('zSeedMin',        5.5);
opts.zThrFoot        = go('zThrFoot',        2.5);
opts.minPeakDist_s   = go('minPeakDist_s',   0.2);
opts.minPeakProm     = go('minPeakProm',     0.75);
opts.minArea         = go('minArea',         2);
opts.eventMaxArea    = go('eventMaxArea',    40);
opts.nmfReps         = go('nmfReps',         5);
opts.nmfMaxIter      = go('nmfMaxIter',      300);
opts.nmfMinPix       = go('nmfMinPix',       4);
opts.r_min           = go('r_min',           3);
opts.r_max           = go('r_max',           8);
opts.r_max_min       = go('r_max_min',       0.5);
opts.r_max_max       = go('r_max_max',       20);
opts.bound           = go('bound',           10);
opts.minEventsSyn    = go('minEventsSyn',    1);
opts.footSigma       = go('footSigma',       0.75);
opts.areaThr         = go('areaThr',         0);
opts.synFilter       = go('synFilter',       []);
opts.saveSynFilter   = go('saveSynFilter',   false);
opts.synFilterFile   = go('synFilterFile',   '');
opts.usePeakWeight   = go('usePeakWeight',   true);
opts.doPlot          = go('doPlot',          true);
opts.uiSeedThresh    = go('uiSeedThresh',    true);
opts.uiClusterSize   = go('uiClusterSize',   true);
opts.uiSynFilter     = go('uiSynFilter',     true);
opts.plotCentroids   = go('plotCentroids',   true);
opts.diagSaveDir     = go('diagSaveDir',     pwd);
opts.diagPrefix      = go('diagPrefix',      'GluDiag');
opts.diagZoom_pad    = go('diagZoom_pad',    15);
opts.diagMaxSyn      = go('diagMaxSyn',      16);
opts.showNMFDiag     = go('showNMFDiag',     true);
opts.nmfDiagMaxEvent = go('nmfDiagMaxEvent', 6);

if ~exist(opts.diagSaveDir,'dir'), mkdir(opts.diagSaveDir); end

%% ===================================================================
%  PROJECT MODE
%% ===================================================================
if isfield(opts,'mode') && strcmpi(opts.mode,'project')
    if nargin < 5 || isempty(S_glu_in)
        error('extractGluSNFR3: project mode requires S_glu_in.');
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
    GluResult = struct('S_glu',S_glu_in,'T_glu',T_glu,...
        'F0_glu',F0_glu,'dFF_glu',dFF_glu,'frameRate',frameRate,'opts',opts);
    return
end

%% ===================================================================
%  STEP 0: PREPROCESS
%% ===================================================================
mov = double(mov2_sub);
[H,W,T] = size(mov);
minPeakDist = max(1, round(opts.minPeakDist_s * frameRate));

mov_sm = imgaussfilt(mov, opts.smoothSigma);

% Robust per-pixel z-score
medImg   = median(mov_sm, 3);
madImg   = mad(mov_sm, 1, 3);
sigmaImg = 1.4826 * madImg + eps;
zMov     = (mov_sm - medImg) ./ sigmaImg;

% Max-projection for UI
zMaxProj = max(zMov, [], 3, 'omitnan');
nzv = zMaxProj(zMaxProj ~= 0);
if ~isempty(nzv), zMaxProj(zMaxProj==0) = median(nzv); end

% Pre-flatten zMov for fast blob-pixel lookup in Step 2
zMovMat = reshape(zMov, H*W, T);   % (H*W) x T

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
    redrawSeed(sld);          % safe: button already created
    uiwait(fig);
    opts.zSeedMin = sld.Value;
    close(fig);
    fprintf('zSeedMin = %.2f\n', opts.zSeedMin);
end

%% ===================================================================
%  STEP 1: PIXEL-PEAK SEEDS
%% ===================================================================
fprintf('Step 1: finding per-pixel peaks...\n');
seedByFrame = cell(T,1);
for y = 1:H
    for x = 1:W
        tr = squeeze(zMov(y,x,:));
        if max(tr) < opts.zSeedMin, continue; end
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
%  STEP 2: BLOB SEGMENTATION + NMF SPLIT OF LARGE BLOBS
%% ===================================================================
% Each detected blob is either kept as one event (small) or factorised
% with NMF(pixels x time) into K spatially distinct sub-events (large).
%
% NMF matrix:  Xmat = blob_pixels x T  (z-score time-traces)
%   Wnmf : blob_pixels x K   — spatial loading per component
%   Hnmf : K x T             — temporal (not stored downstream)
% Each pixel → argmax(Wnmf, dim=2) → component assignment.
%
% The resulting eventList (possibly larger than blob count) feeds Step 3.
fprintf('Step 2: segmenting blobs (NMF-split threshold: %d px)...\n', opts.eventMaxArea);

eventList    = struct('frame',{},'pixIdx',{},'centroid_xy',{},...
                      'peakZ',{},'pixVal',{},'meanZ',{},'nmfSplit',{});
compCount    = 0;
nmfDiagData  = struct('frame',{},'origPix',{},'K',{},'Wnmf',{},'subCents',{});
nmfDiagCount = 0;

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

        % --- Small blob: single event ---
        if numel(pixIdx) <= opts.eventMaxArea
            compCount = compCount + 1;
            eventList(compCount) = buildEvent(t,pixIdx,c,r,pixVal,false);
            continue
        end

        % --- Large blob: NMF split ---
        K    = max(2, round(numel(pixIdx)/opts.eventMaxArea));
        Xmat = max(zMovMat(double(pixIdx(:)),:), 0);   % blob_px x T

        usedNMF = false;
        Wnmf    = [];
        if size(Xmat,1)>=K && size(Xmat,2)>=K
            try
                  optnmf = statset('MaxIter',5);
                [Wnmf,~] = nnmf(Xmat, K,...
                    'replicates',opts.nmfReps,...
                    'Options',   optnmf);
                usedNMF = true;
            catch
            end
        end

        if ~usedNMF
            compCount = compCount + 1;
            eventList(compCount) = buildEvent(t,pixIdx,c,r,pixVal,false);
            continue
        end

        [~,pixComp] = max(Wnmf,[],2);   % blob_px x 1, values 1..K
        subCents = nan(K,2);
        nAdded   = 0;

        for k = 1:K
            mask   = (pixComp==k);
            subPix = pixIdx(mask);
            if numel(subPix) < opts.nmfMinPix, continue; end
            [sr,sc_] = ind2sub([H,W], subPix);
            subVal   = pixVal(mask);
            subCents(k,:) = [mean(sc_) mean(sr)];
            compCount = compCount + 1;
            eventList(compCount) = buildEvent(t,subPix,sc_,sr,subVal,true);
            nAdded = nAdded + 1;
        end

        if nAdded == 0
            % NMF produced nothing large enough: keep original blob
            compCount = compCount + 1;
            eventList(compCount) = buildEvent(t,pixIdx,c,r,pixVal,false);
        end

        if nmfDiagCount < opts.nmfDiagMaxEvent
            nmfDiagCount = nmfDiagCount + 1;
            nmfDiagData(nmfDiagCount).frame    = t;
            nmfDiagData(nmfDiagCount).origPix  = pixIdx;
            nmfDiagData(nmfDiagCount).K        = K;
            nmfDiagData(nmfDiagCount).Wnmf     = Wnmf;
            nmfDiagData(nmfDiagCount).subCents = subCents;
        end
    end
end

nSplit = sum([eventList.nmfSplit]);
fprintf('  %d events extracted (%d from NMF splits).\n', numel(eventList), nSplit);

% NMF blob-split gallery — one figure, each row = one split event
if opts.doPlot && opts.showNMFDiag && nmfDiagCount > 0
    pad = opts.diagZoom_pad;

    % Max K across all splits (sets number of columns)
    maxK = max([nmfDiagData(1:nmfDiagCount).K]);
    nCols_gal = maxK + 2;   % original blob | K components | centroid map
    nRows_gal = nmfDiagCount;

    fh_gal = figure('Name','NMF blob-split gallery',...
                    'Position',[40 40 260*nCols_gal  250*nRows_gal]);
    tiledlayout(nRows_gal, nCols_gal, ...
        'TileSpacing','compact','Padding','compact');

    for di = 1:nmfDiagCount
        d = nmfDiagData(di);
        [pr,pc_] = ind2sub([H,W], d.origPix);
        xc_ = round(mean(pc_)); yc_ = round(mean(pr));
        xrng = max(1,xc_-pad):min(W,xc_+pad);
        yrng = max(1,yc_-pad):min(H,yc_+pad);
        cmap_k = lines(d.K);

        blobImg = zeros(H,W);
        blobImg(d.origPix) = zMov(d.origPix + (d.frame-1)*H*W);

        % Column 1: original blob at event frame
        ax0 = nexttile;
        imagesc(ax0, blobImg(yrng,xrng));
        axis(ax0,'image'); axis(ax0,'off');
        colormap(ax0, turbo);
        if di == 1
            title(ax0, 'Original blob', 'FontSize',8);
        end
        ylabel(ax0, sprintf('Split %d\nt=%d K=%d', di, d.frame, d.K), ...
               'FontSize',7, 'Rotation',0, 'HorizontalAlignment','right');

        % Columns 2..K+1: NMF spatial components
        for k = 1:maxK
            ax_k = nexttile;
            if k <= d.K
                Wk = zeros(H,W);
                Wk(d.origPix) = d.Wnmf(:,k);
                imagesc(ax_k, Wk(yrng,xrng));
                axis(ax_k,'image'); axis(ax_k,'off');
                colormap(ax_k, hot);
                rectangle(ax_k,'Position',[0.5 0.5 numel(xrng) numel(yrng)],...
                    'EdgeColor',cmap_k(k,:),'LineWidth',2);
                if di == 1
                    title(ax_k, sprintf('Comp %d',k), 'FontSize',8);
                end
            else
                % Empty tile for splits with fewer components than maxK
                axis(ax_k,'off');
            end
        end

        % Last column: sub-event centroids overlaid on blob
        ax_s = nexttile;
        imagesc(ax_s, blobImg(yrng,xrng));
        axis(ax_s,'image'); axis(ax_s,'off');
        colormap(ax_s, gray);
        hold(ax_s,'on');
        for k = 1:d.K
            sc_ = d.subCents(k,:);
            if any(isnan(sc_)), continue; end
            plot(ax_s, sc_(1)-xrng(1)+1, sc_(2)-yrng(1)+1, ...
                '+', 'Color',cmap_k(k,:), 'MarkerSize',12, 'LineWidth',2);
        end
        hold(ax_s,'off');
        if di == 1
            title(ax_s, 'Sub-centroids', 'FontSize',8);
        end
    end

    sgtitle(sprintf('NMF blob-split gallery  (%d splits, threshold=%d px)', ...
            nmfDiagCount, opts.eventMaxArea), 'FontSize',11);
    saveDiag(fh_gal, opts, 'NMF_split_gallery');
end

if isempty(eventList)
    warning('No events detected. Returning empty result.');
    GluResult = struct('eventList',[],'S_glu',[],'T_glu',[],...
                       'dFF_glu',[],'opts',opts,'frameRate',frameRate);
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
    fig = uifigure('Name','Adjust r\_max — max synapse diameter before k-means split',...
                   'Position',[80 80 920 590]);
    ax  = uiaxes(fig,'Position',[20 115 880 445]);
    axis(ax,'ij'); axis(ax,'equal'); xlim(ax,[1 W]); ylim(ax,[1 H]);
    title(ax,'Event centroids (circles=normal, diamonds=NMF sub-events)');
    hold(ax,'on');
    b = opts.bound;
    rectangle(ax,'Position',[b+1,b+1,W-2*b,H-2*b],...
        'EdgeColor','w','LineWidth',2,'LineStyle','--');
    sld = uislider(fig,'Position',[20 78 880 3],...
                   'Limits',[opts.r_max_min opts.r_max_max],'Value',opts.r_max);
    lbl = uilabel(fig,'Position',[20 50 600 26],...
                  'Text',sprintf('r_max = %.2f px',opts.r_max),'FontSize',13);
    btn = uibutton(fig,'Text','Accept','Position',[720 14 180 36],...
                   'FontSize',13,'ButtonPushedFcn',@(~,~) uiresume(fig)); %#ok<NASGU>
    fig.UserData.ax = ax;
    fig.UserData.sh = [];
    sld.ValueChangedFcn = @(src,~) redrawCluster(fig,lbl,src.Value);
    redrawCluster(fig,lbl,sld.Value);   % safe: button already exists
    uiwait(fig);
    opts.r_max = sld.Value;
    close(fig);
    fprintf('r_max = %.2f px\n', opts.r_max);
end

% r_min graph connectivity
r2 = opts.r_min^2;
A  = false(Nevent,Nevent);
for i = 1:Nevent
    dx = cent(i,1)-cent(:,1); dy = cent(i,2)-cent(:,2);
    A(i,:) = (dx.^2+dy.^2) <= r2;
end
clusters = accumarray(conncomp(graph(A))',(1:Nevent)',[],@(v){v});

synID      = zeros(Nevent,1);
synCounter = 0;

for c = 1:numel(clusters)
    idx = clusters{c};
    pts = cent(idx,:);
    bboxDiag = hypot(max(pts(:,1))-min(pts(:,1)), max(pts(:,2))-min(pts(:,2)));

    if bboxDiag <= opts.r_max || numel(idx) <= 2
        synCounter = synCounter + 1;
        synID(idx) = synCounter;
        continue
    end

    K = max(2, ceil(bboxDiag/opts.r_max));
    try
        lab = kmeans(pts,K,'Replicates',5,'MaxIter',300);
    catch
        synCounter = synCounter + 1;
        synID(idx) = synCounter;
        continue
    end
    for k = 1:K
        si = idx(lab==k);
        if isempty(si), continue; end
        synCounter = synCounter + 1;
        synID(si) = synCounter;
    end
end

fprintf('  %d events → %d synapses (r_min=%.1f, r_max=%.1f px)\n',...
    Nevent, synCounter, opts.r_min, opts.r_max);

if opts.doPlot && opts.plotCentroids && synCounter > 0
    fh = figure('Name','Event centroids by synapse');
    hold on; axis ij; axis equal; xlim([1 W]); ylim([1 H]);
    b = opts.bound;
    rectangle('Position',[b+1,b+1,W-2*b,H-2*b],...
        'EdgeColor','w','LineWidth',2,'LineStyle','--');
    cmap_s = hsv(synCounter);
    for s = 1:synCounter
        ei = find(synID==s); if isempty(ei), continue; end
        isSplit = [eventList(ei).nmfSplit];
        if any(~isSplit)
            scatter(cent(ei(~isSplit),1),cent(ei(~isSplit),2),...
                14,cmap_s(s,:),'filled');
        end
        if any(isSplit)
            scatter(cent(ei(isSplit),1),cent(ei(isSplit),2),...
                22,cmap_s(s,:),'d','filled');
        end
    end
    title('Event centroids  (circles = normal, diamonds = NMF sub-events)');
    set(gca,'Color','k'); hold off;
    saveDiag(fh,opts,'centroids_by_synapse');
end

%% ===================================================================
%  STEP 4: BUILD S_glu — mean accumulated z-score footprints
%% ===================================================================
% For each event, only pixels within r_max of the synapse centroid are
% used. This is critical for NMF sub-events, whose assigned pixels can
% be scattered across the original large blob — cropping to r_max gives
% each sub-synapse a clean, localised footprint.
% Normal events are also cropped, which removes any stray pixels at the
% edge of a blob that belong spatially to a neighbouring synapse.
fprintf('Step 4: building synapse footprints (r_max crop = %.1f px)...\n', opts.r_max);

% Pre-compute synapse centroids from clustering (weighted mean of event centroids)
% We need these before building footprints so we can crop relative to them.
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

    % --- Crop pixels to within r_max of this synapse's centroid ---
    cx = synCentroid0(s,1);   % synapse centroid x (col)
    cy = synCentroid0(s,2);   % synapse centroid y (row)
    [pr, pc] = ind2sub([H,W], pixIdx_d);
    inRadius = (pc - cx).^2 + (pr - cy).^2 <= opts.r_max^2;
    pixIdx_d = pixIdx_d(inRadius);
    pixVal_d = pixVal_d(inRadius);

    if isempty(pixIdx_d), continue; end   % entire event outside r_max (rare)

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

% Footprint features
synArea = zeros(newNsyn,1);
synAmp  = zeros(newNsyn,1);
for s = 1:newNsyn
    fp = S_glu(:,:,s);
    nz = fp(fp > opts.areaThr);
    synArea(s) = numel(nz);
    synAmp(s)  = ifelse(isempty(nz), 0, mean(nz));
end

% Quick amplitude proxy: dot zMov onto each footprint
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

% Interactive filter UI
if opts.doPlot && opts.uiSynFilter && newNsyn > 0
    fig = uifigure('Name','Synapse filter','Position',[80 80 940 590]);
    ax  = uiaxes(fig,'Position',[20 235 900 330]);
    axis(ax,'ij'); axis(ax,'equal'); xlim(ax,[1 W]); ylim(ax,[1 H]);
    title(ax,'Red = kept     Green = filtered');
    imagesc(ax, zMaxProj); colormap(ax,gray(256)); hold(ax,'on');
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
    redrawFilter();   % safe: button already exists
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

duration_s      = T / frameRate;
synEventRate_Hz = synEventCount / max(duration_s,eps);
fprintf('  Kept %d synapses.  Median rate: %.3f Hz\n', newNsyn, median(synEventRate_Hz));

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
    % Event count map
    fh = figure('Name','Event count map');
    imagesc(Smap_counts); axis image; axis off; colormap hot; colorbar;
    title('Event count map — kept synapses');
    saveDiag(fh,opts,'event_count_map');

    % Footprint overview — show_footprnt from your library
    fh = figure('Name','Synapse footprints');
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    show_footprnt(S_glu, mean(mov2_sub,3));
    sgtitle('Synapse footprints (turbo) | overlay on mean movie');
    saveDiag(fh,opts,'synapse_footprints');

    % Event-rate heatmap — show_footprint_heatmap from your library
    fh = figure('Name','Footprints colored by event rate');
    tiledlayout(1,1);
    show_footprint_heatmap(S_glu, synEventRate_Hz);
    title('Footprints colored by event rate (Hz)');
    saveDiag(fh,opts,'synapse_footprints_eventrate');

    % Build RGB composite for output field
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
%  STEP 7: COMPUTE T_glu  (least-squares projection)
%% ===================================================================
fprintf('Step 7: computing T_glu...\n');
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
    tvec  = (0:T-1)/frameRate;
    pad   = opts.diagZoom_pad;
    nShow = min(newNsyn, opts.diagMaxSyn);
    nCols = 4; nRows = ceil(nShow/2);

    % (A) Per-synapse: zoomed footprint + trace
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
        xlabel(ax2,'Time (s)','FontSize',7);
        title(ax2,sprintf('n=%d  %.2f Hz',...
            synEventCount(s),synEventRate_Hz(s)),'FontSize',8);
        set(ax2,'FontSize',7);
    end
    sgtitle('Per-synapse footprint (left) | T\_glu trace (right)','FontSize',11);
    saveDiag(fh,opts,'per_synapse_footprint_trace');

    % (C) All footprints grid
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
GluResult = struct();
GluResult.eventList       = eventList;
GluResult.Smap_event      = Smap_counts;
GluResult.S_glu           = S_glu;
GluResult.synCentroid     = synCentroid;
GluResult.synEventCount   = synEventCount;
GluResult.synEventRate_Hz = synEventRate_Hz;
GluResult.T_glu           = T_glu;
GluResult.F0_glu          = F0_glu;
GluResult.dFF_glu         = dFF_glu;
GluResult.RGB             = RGB;
GluResult.opts            = opts;
GluResult.frameRate       = frameRate;

%% ===================================================================
%  LOCAL HELPERS  (nested functions — share workspace via closure)
%% ===================================================================

    % Build one eventList entry
    function ev = buildEvent(t, pixIdx, c, r, pixVal, isSplit)
        ev.frame       = t;
        ev.pixIdx      = pixIdx;
        ev.centroid_xy = [mean(c) mean(r)];
        ev.peakZ       = max(pixVal);
        ev.pixVal      = pixVal;
        ev.meanZ       = mean(pixVal);
        ev.nmfSplit    = isSplit;
    end

    % Save figure to diagSaveDir as PNG
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

    % Ternary helper
    function v = ifelse(cond, a, b)
        if cond; v=a; else; v=b; end
    end

    % Percentile-based slider limits
    function L = prcLims(v, pLo, pHi, epsPad)
        v = v(isfinite(v(:)));
        if isempty(v), L=[0 1]; return; end
        L = [prctile(v,pLo) prctile(v,pHi)];
        if L(1)==L(2), L=[L(1) L(1)+epsPad]; end
    end

    % Seed threshold preview
    function redrawSeed(src)
        zthr = src.Value;
        lbl.Text = sprintf('zSeedMin = %.2f', zthr);
        imagesc(ax1, max(zMov>zthr,[],3)); colormap(ax1,gray);
        axis(ax1,'image'); axis(ax1,'off');
        title(ax1,sprintf('max(zMov > %.2f, [], 3)',zthr),'FontSize',9);
        imagesc(ax2, zMaxProj); colormap(ax2,turbo);
        axis(ax2,'image'); axis(ax2,'off');
        title(ax2,'Max z-score projection','FontSize',9);
        drawnow;
    end

    % r_max cluster preview — recomputes full clustering with current rmax
    function redrawCluster(figH, lblH, rmax)
        lblH.Text = sprintf('r_max = %.2f px', rmax);
        centL    = cent;
        N_       = size(centL,1);
        r2L      = opts.r_min^2;
        axL      = figH.UserData.ax;
        isSplit_ = [eventList.nmfSplit];

        % r_min graph connectivity
        AL = false(N_);
        for ii = 1:N_
            dx_ = centL(ii,1)-centL(:,1);
            dy_ = centL(ii,2)-centL(:,2);
            AL(ii,:) = (dx_.^2+dy_.^2) <= r2L;
        end
        cls_ = accumarray(conncomp(graph(AL))',(1:N_)',[],@(v){v});

        % Split large clusters with k-means (same logic as Step 3)
        synID_ = zeros(N_,1);
        sc_    = 0;
        for ci_ = 1:numel(cls_)
            ii_  = cls_{ci_};
            pts_ = centL(ii_,:);
            bd_  = hypot(max(pts_(:,1))-min(pts_(:,1)),...
                         max(pts_(:,2))-min(pts_(:,2)));
            if bd_ <= rmax || numel(ii_) <= 2
                sc_ = sc_+1; synID_(ii_) = sc_;
                continue
            end
            K__ = max(2, ceil(bd_/rmax));
            try
                lb_ = kmeans(pts_,K__,'Replicates',3,'MaxIter',200);
            catch
                sc_ = sc_+1; synID_(ii_) = sc_;
                continue
            end
            for kk_ = 1:K__
                si_ = ii_(lb_==kk_);
                if isempty(si_), continue; end
                sc_ = sc_+1; synID_(si_) = sc_;
            end
        end

        % Assign colors by cluster ID
        cm  = hsv(max(sc_,1));
        % Shuffle so adjacent IDs don't share similar colors
        rng(0); perm = randperm(max(sc_,1));
        cm  = cm(perm,:);
        C   = cm(synID_,:);

        % Update or create scatter handles (circles + diamonds)
        if isempty(figH.UserData.sh) || ...
                ~all(isgraphics(figH.UserData.sh))
            cla(axL);   % clear previous points
            hold(axL,'on');
            h1 = scatter(axL, centL(~isSplit_,1), centL(~isSplit_,2), ...
                         14, C(~isSplit_,:), 'filled', 'o');
            h2 = scatter(axL, centL(isSplit_,1),  centL(isSplit_,2), ...
                         22, C(isSplit_,:),  'filled', 'd');
            hold(axL,'off');
            figH.UserData.sh = [h1 h2];
        else
            figH.UserData.sh(1).XData = centL(~isSplit_,1);
            figH.UserData.sh(1).YData = centL(~isSplit_,2);
            figH.UserData.sh(1).CData = C(~isSplit_,:);
            if ~isempty(centL(isSplit_,1))
                figH.UserData.sh(2).XData = centL(isSplit_,1);
                figH.UserData.sh(2).YData = centL(isSplit_,2);
                figH.UserData.sh(2).CData = C(isSplit_,:);
            end
        end
        title(axL, sprintf('r_{max}=%.1f px | %d synapses', rmax, sc_), ...
              'FontSize',10);
        drawnow;
    end

    % Synapse filter preview
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

    % Min/max slider pair
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

end % extractGluSNFR3