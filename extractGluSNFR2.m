function GluResult = extractGluSNFR2(mov2_sub, mov2_mc_filt, exposuretime2, opts, S_glu_in)
% extractGluSNFR  (MODIFIED)
%   Detect pixel-wise iGluSNFR events, cluster into synapses, build S_glu (H×W×Nsyn),
%   remove boundary synapses, visualize, and compute T_glu by least-squares projection.
%
% MODIFICATIONS vs original:
%   1) S_glu built from accumulated z-score footprint masks (not centroid impulses).
%   2) Large single-frame event blobs are split by NMF BEFORE centroid clustering.
%      If a segmented event region exceeds opts.eventMaxArea pixels, NMF is run on
%      that blob's pixels × time (zMov) to find K spatial sub-sources. Each
%      sub-source becomes its own entry in eventList with its own centroid.
%      These expanded centroids then feed the normal r_min/r_max clustering.
%   3) Centroid clustering restored to original r_min graph + k-means (no NMF there).
%   4) Comprehensive diagnostic plots saved to disk and shown interactively.
%
% USAGE
%   GluResult = extractGluSNFR(mov2_sub, mov2_mc_filt, exposuretime2);
%   GluResult = extractGluSNFR(mov2_sub, mov2_mc_filt, exposuretime2, opts);
%
% KEY opts fields (new / changed)
%   opts.diagSaveDir   : folder to save diagnostic PNGs (default: pwd)
%   opts.diagPrefix    : filename prefix for saved figures (default: 'GluDiag')
%   opts.nmfMaxIter    : NMF iterations for large-cluster splitting (default: 200)
%   opts.nmfReps       : NMF replicates (default: 5)
%   opts.splitBySpatialNMF : if true (default), use spatial NMF on event masks
%                            to split clusters; otherwise fall back to centroid k-means
%   opts.diagZoom_pad  : pixels around synapse centroid for zoomed footprint (default: 15)
%   opts.diagMaxSyn    : max synapses shown in the per-synapse diagnostic (default: 16)

%% ========================= DEFAULTS =========================
if nargin < 4 || isempty(opts), opts = struct; end
    function val = getOpt(o, field, default)
        if isfield(o, field) && ~isempty(o.(field))
            val = o.(field);
        else
            val = default;
        end
    end

frameRate = 1/exposuretime2;

% z-score / peak seed params
opts.smoothSigma      = getOpt(opts,'smoothSigma',      1);
opts.zSeedMin         = getOpt(opts,'zSeedMin',         3);
opts.zThrFoot         = getOpt(opts,'zThrFoot',         2.5);
opts.minPeakDist_s    = getOpt(opts,'minPeakDist_s',    0.2);
opts.minPeakProm      = getOpt(opts,'minPeakProm',      0.75);
opts.minArea          = getOpt(opts,'minArea',          2);
opts.eventPixMax      = getOpt(opts,'eventPixMax',      5);
opts.r_min            = getOpt(opts,'r_min',            3);    % px: graph connectivity radius for clustering
opts.r_max            = getOpt(opts,'r_max',            8);    % px: max cluster size before k-means split
opts.bound            = getOpt(opts,'bound',            10);
opts.doPlot           = getOpt(opts,'doPlot',           true);
opts.uiSeedThresh     = getOpt(opts,'uiSeedThresh',     true);
opts.uiSynFilter      = getOpt(opts,'uiSynFilter',      true);
opts.plotCentroids    = getOpt(opts,'plotCentroids',    true);
opts.usePeakWeight    = getOpt(opts,'usePeakWeight',    true);
opts.minEventsSyn     = getOpt(opts,'minEventsSyn',     3);
opts.uiClusterSize    = getOpt(opts,'uiClusterSize',    true);
opts.r_max_min        = getOpt(opts,'r_max_min',        0.001);
opts.r_max_max        = getOpt(opts,'r_max_max',        10);
opts.footSigma        = getOpt(opts,'footSigma',        1);
opts.areaThr          = getOpt(opts,'areaThr',          0);
opts.synFilter        = getOpt(opts,'synFilter',        []);
opts.saveSynFilter    = getOpt(opts,'saveSynFilter',    false);
opts.synFilterFile    = getOpt(opts,'synFilterFile',    '');

% --- NMF per-event blob splitting (Step 2, before clustering) ---
opts.diagSaveDir      = getOpt(opts,'diagSaveDir',      pwd);
opts.diagPrefix       = getOpt(opts,'diagPrefix',       'GluDiag');
opts.eventMaxArea     = getOpt(opts,'eventMaxArea',     40);   % px^2: blobs larger than this get NMF-split
opts.nmfMaxIter       = getOpt(opts,'nmfMaxIter',       200);
opts.nmfReps          = getOpt(opts,'nmfReps',          5);
opts.nmfMinPix        = getOpt(opts,'nmfMinPix',        4);    % min pixels per NMF sub-event to keep
opts.diagZoom_pad     = getOpt(opts,'diagZoom_pad',     15);
opts.diagMaxSyn       = getOpt(opts,'diagMaxSyn',       16);
opts.showNMFDiag      = getOpt(opts,'showNMFDiag',      true);
opts.nmfDiagMaxEvent  = getOpt(opts,'nmfDiagMaxEvent',  6);    % max split events to show diagnostic for

if ~exist(opts.diagSaveDir,'dir')
    mkdir(opts.diagSaveDir);
end

%%  ========================= MODE SWITCH =========================
if ~isfield(opts,'mode') || isempty(opts.mode)
    opts.mode = 'detect';
end

if strcmpi(opts.mode, 'project')
    if nargin < 5 || isempty(S_glu_in)
        error('Project mode requires S_glu_in (H x W x Nsyn).');
    end
    mov = double(mov2_sub);
    [H,W,T] = size(mov);
    if ~isempty(mov2_mc_filt)
        sz1 = size(mov2_sub); sz2 = size(mov2_mc_filt);
        if any(sz1 ~= sz2)
            error('Size mismatch between mov2_sub and mov2_mc_filt.');
        end
    end
    [Hs,Ws,Nsyn] = size(S_glu_in);
    assert(H==Hs && W==Ws, 'S_glu_in must match movie spatial size.');
    Npix = H*W;
    S = reshape(double(S_glu_in), Npix, Nsyn);
    M = reshape(mov, Npix, T);
    M(~isfinite(M)) = 0;
    T_glu = (S \ M);
    F0_glu = []; dFF_glu = [];
    if ~isempty(mov2_mc_filt)
        M0 = reshape(double(mov2_mc_filt), Npix, T);
        M0(~isfinite(M0)) = 0;
        F0_glu  = (S \ M0);
        dFF_glu = T_glu ./ mean(F0_glu, 2);
    end
    GluResult = struct('S_glu',S_glu_in,'T_glu',T_glu,...
        'F0_glu',F0_glu,'dFF_glu',dFF_glu,...
        'frameRate',1/exposuretime2,'opts',opts);
    return
end

%% ========================= PREPROCESSING =========================
mov = double(mov2_sub);
[H,W,T] = size(mov);
minPeakDist = max(1, round(opts.minPeakDist_s * frameRate));
mov = imgaussfilt(mov, opts.smoothSigma);

% Robust per-pixel z-score
medImg   = median(mov, 3);
madImg   = mad(mov, 1, 3);
sigmaImg = 1.4826 * madImg + eps;
zMov     = (mov - medImg) ./ sigmaImg;

zMaxProj = max(zMov, [], 3, 'omitnan');
zMaxProj(zMaxProj==0) = median(tovec(zMaxProj(zMaxProj~=0)));

%% ========================= UI: SEED THRESHOLD =========================
if opts.doPlot && opts.uiSeedThresh
    zSeedMin0 = opts.zSeedMin;
    zMin = 0;
    zMax = max(6, ceil(max(zMov(:),[],'omitnan')));

    % Layout (all positions in pixels, [left bottom width height]):
    %   Figure: 900 wide x 560 tall
    %   Axes:   top portion, 420 tall each, starting at y=120
    %   Slider: y=75 (needs ~30px height for hit area)
    %   Label:  y=50
    %   Button: y=15
    fig = uifigure('Name','Adjust zSeedMin','Position',[100 100 900 560]);
    ax  = uiaxes(fig,'Position',[20  120  420  400]);
    ax2 = uiaxes(fig,'Position',[460 120  420  400]);
    sld = uislider(fig,'Position',[20 80 860 3], ...
                   'Limits',[zMin zMax],'Value',zSeedMin0);
    lbl = uilabel(fig,'Position',[20 50 500 22], ...
                  'Text',sprintf('zSeedMin = %.2f', sld.Value));
    btn = uibutton(fig,'Text','Accept zSeedMin', ...
                   'Position',[700 12 180 32], ...
                   'ButtonPushedFcn', @(~,~) uiresume(fig));

    % Wire callback and do initial draw AFTER button exists
    sld.ValueChangedFcn = @(src,~) redraw_slider(src);
    redraw_slider(sld);

    uiwait(fig);
    opts.zSeedMin = sld.Value;
    close(fig);
    fprintf('zSeedMin updated to %.2f\n', opts.zSeedMin);
end

%% ========================= STEP 1: PIXEL PEAK SEEDS =========================
fprintf('Finding per-pixel peaks...\n');
seedByFrame = cell(T,1);
for y = 1:H
    for x = 1:W
        tr = squeeze(zMov(y,x,:));
        if max(tr) < opts.zSeedMin, continue; end
        [pks, locs] = findpeaks(tr, ...
            'MinPeakDistance',   minPeakDist, ...
            'MinPeakProminence', opts.minPeakProm);
        if isempty(locs), continue; end
        locs = locs(pks >= opts.zSeedMin);
        if isempty(locs), continue; end
        lin = sub2ind([H,W], y, x);
        for k = 1:numel(locs)
            seedByFrame{locs(k)}(end+1) = lin; %#ok<AGROW>
        end
    end
end
fprintf('Total seed peaks: %d\n', sum(cellfun(@numel, seedByFrame)));

%% ========================= STEP 2: SEGMENT EVENT FOOTPRINTS (+ NMF split large blobs) =========================
% For each segmented event blob:
%   - If numel(pixIdx) <= opts.eventMaxArea  →  add as single event (original behaviour)
%   - If numel(pixIdx) >  opts.eventMaxArea  →  run NMF on pixels × time (zMov) to
%     separate spatially distinct sub-sources; each component becomes its own event.
% The resulting eventList (possibly expanded) feeds the clustering step unchanged.

eventList = struct('frame',{},'pixIdx',{},'centroid_xy',{},...
                   'peakZ',{},'pixVal',{},'meanZ',{},'nmfSplit',{});
compCount    = 0;
nmfDiagData  = struct('frame',{},'origPixIdx',{},'K',{},'Wnmf',{},'subCents',{});
nmfDiagCount = 0;

for t = 1:T
    seeds = seedByFrame{t};
    if isempty(seeds), continue; end
    zFrame = zMov(:,:,t);
    bw = zFrame > opts.zThrFoot;
    bw = bwareaopen(bw, opts.minArea);
    bw = imopen(bw, strel('disk',1));
    bw = imfill(bw,'holes');
    CC = bwconncomp(bw, 8);
    if CC.NumObjects == 0, continue; end
    L = labelmatrix(CC);
    seedLabels = unique(L(seeds));
    seedLabels(seedLabels==0) = [];
    if isempty(seedLabels), continue; end

    for lab = seedLabels(:)'
        pixIdx = find(L == lab);
        if numel(pixIdx) < opts.minArea, continue; end
        [r,c]  = ind2sub([H,W], pixIdx);
        pixVal = zFrame(pixIdx);

        % ---- Small blob: single event (original behaviour) ----
        if numel(pixIdx) <= opts.eventMaxArea
            compCount = compCount + 1;
            eventList(compCount).frame        = t;
            eventList(compCount).pixIdx       = pixIdx;
            eventList(compCount).centroid_xy  = [mean(c) mean(r)];
            eventList(compCount).peakZ        = max(pixVal);
            eventList(compCount).pixVal       = pixVal;
            eventList(compCount).meanZ        = mean(pixVal);
            eventList(compCount).nmfSplit     = false;
            continue
        end

        % ---- Large blob: NMF split into K sub-events ----
        % Xmat = pixels-in-blob × T  (z-score time traces for each blob pixel)
        % NMF: Xmat ≈ Wnmf (pixels×K) * Hnmf (K×T)
        % Each column of Wnmf is a spatial component (sub-synapse candidate).
        % We threshold each component to get sub-event pixIdx and centroid.

        % Estimate K from blob area
        K = max(2, round(numel(pixIdx) / opts.eventMaxArea));

        % Build pixel × time matrix for this blob
        Xmat = max(squeeze(zMov(pixIdx + 0*H*W)), 0);  % numel(pixIdx) × 1, wrong — need all T
        % Correct: index each frame for the blob pixels
        zMovMat = reshape(zMov, H*W, T);               % (H*W) × T
        Xmat    = max(zMovMat(double(pixIdx(:)), :), 0); % Npix_blob × T

        % Run NMF
        usedNMF = false;
        Wnmf    = [];
        if size(Xmat,1) >= K && size(Xmat,2) >= K
            try
                [Wnmf, ~] = nnmf(Xmat, K, ...
                    'replicates', opts.nmfReps, ...
                    'maxiter',    opts.nmfMaxIter);
                usedNMF = true;
            catch
                usedNMF = false;
            end
        end

        if ~usedNMF
            % Fallback: treat as single event
            compCount = compCount + 1;
            eventList(compCount).frame        = t;
            eventList(compCount).pixIdx       = pixIdx;
            eventList(compCount).centroid_xy  = [mean(c) mean(r)];
            eventList(compCount).peakZ        = max(pixVal);
            eventList(compCount).pixVal       = pixVal;
            eventList(compCount).meanZ        = mean(pixVal);
            eventList(compCount).nmfSplit     = false;
            continue
        end

        % Assign each pixel to the NMF component with highest spatial loading
        [~, pixComp] = max(Wnmf, [], 2);   % Npix_blob × 1, values 1..K

        subCents = nan(K, 2);
        nAdded   = 0;
        for k = 1:K
            subMask = pixComp == k;
            subPix  = pixIdx(subMask);
            if numel(subPix) < opts.nmfMinPix, continue; end
            [sr,sc_] = ind2sub([H,W], subPix);
            subVal   = pixVal(subMask);
            subCents(k,:) = [mean(sc_) mean(sr)];

            compCount = compCount + 1;
            eventList(compCount).frame        = t;
            eventList(compCount).pixIdx       = subPix;
            eventList(compCount).centroid_xy  = [mean(sc_) mean(sr)];
            eventList(compCount).peakZ        = max(subVal);
            eventList(compCount).pixVal       = subVal;
            eventList(compCount).meanZ        = mean(subVal);
            eventList(compCount).nmfSplit     = true;
            nAdded = nAdded + 1;
        end

        % If NMF produced nothing useful, add original blob as single event
        if nAdded == 0
            compCount = compCount + 1;
            eventList(compCount).frame        = t;
            eventList(compCount).pixIdx       = pixIdx;
            eventList(compCount).centroid_xy  = [mean(c) mean(r)];
            eventList(compCount).peakZ        = max(pixVal);
            eventList(compCount).pixVal       = pixVal;
            eventList(compCount).meanZ        = mean(pixVal);
            eventList(compCount).nmfSplit     = false;
        end

        % Store diagnostic info for the first nmfDiagMaxEvent splits
        if nmfDiagCount < opts.nmfDiagMaxEvent
            nmfDiagCount = nmfDiagCount + 1;
            nmfDiagData(nmfDiagCount).frame      = t;
            nmfDiagData(nmfDiagCount).origPixIdx = pixIdx;
            nmfDiagData(nmfDiagCount).K          = K;
            nmfDiagData(nmfDiagCount).Wnmf       = Wnmf;
            nmfDiagData(nmfDiagCount).subCents   = subCents;
        end
    end
end

nSplit = sum([eventList.nmfSplit]);
fprintf('Extracted %d event regions (%d from NMF blob splits).\n', numel(eventList), nSplit);

% ---- NMF blob-split diagnostic figures ----
if opts.doPlot && opts.showNMFDiag && nmfDiagCount > 0
    fprintf('Generating NMF blob-split diagnostic for %d events...\n', nmfDiagCount);
    pad = opts.diagZoom_pad;
    for di = 1:nmfDiagCount
        d    = nmfDiagData(di);
        K    = d.K;
        pix  = d.origPixIdx;
        [pr,pc] = ind2sub([H,W], pix);
        xc   = round(mean(pc)); yc = round(mean(pr));
        xrng = max(1,xc-pad):min(W,xc+pad);
        yrng = max(1,yc-pad):min(H,yc+pad);

        fh_nmf = figure('Name', sprintf('NMF blob split — event %d (t=%d, K=%d)', di, d.frame, K), ...
                        'Position', [60 60 300*(K+2) 320]);
        cmap_k = lines(K);

        % Panel 1: original blob z-score at event frame
        ax0 = subplot(1, K+2, 1);
        zFr = zMov(:,:,d.frame);
        blobImg = zeros(H,W);
        blobImg(pix) = zFr(pix);
        imagesc(ax0, blobImg(yrng,xrng)); axis(ax0,'image'); axis(ax0,'off');
        colormap(ax0, turbo); title(ax0, sprintf('Original blob\n(t=%d)',d.frame), 'FontSize',8);

        % Panels 2..K+1: NMF spatial components
        for k = 1:K
            Wslice = zeros(H,W);
            Wslice(pix) = d.Wnmf(:,k);
            ax_k = subplot(1, K+2, k+1);
            imagesc(ax_k, Wslice(yrng,xrng)); axis(ax_k,'image'); axis(ax_k,'off');
            colormap(ax_k, hot);
            title(ax_k, sprintf('NMF comp %d', k), 'FontSize',8);
            rectangle(ax_k, 'Position',[0.5 0.5 numel(xrng) numel(yrng)], ...
                'EdgeColor', cmap_k(k,:), 'LineWidth', 2);
        end

        % Last panel: sub-event centroids overlaid on blob
        ax_s = subplot(1, K+2, K+2);
        imagesc(ax_s, blobImg(yrng,xrng)); axis(ax_s,'image'); axis(ax_s,'off');
        colormap(ax_s, gray); hold(ax_s,'on');
        for k = 1:K
            sc_ = d.subCents(k,:);
            if any(isnan(sc_)), continue; end
            plot(ax_s, sc_(1)-xrng(1)+1, sc_(2)-yrng(1)+1, '+', ...
                'Color', cmap_k(k,:), 'MarkerSize',10, 'LineWidth',2);
        end
        hold(ax_s,'off');
        title(ax_s, 'Sub-centroids', 'FontSize',8);

        sgtitle(sprintf('NMF blob split %d | frame %d | K=%d components', di, d.frame, K), 'FontSize',10);
        saveDiag(fh_nmf, opts, sprintf('NMF_blobsplit_%03d', di));
    end
end

if isempty(eventList)
    warning('No events detected. Returning empty results.');
    GluResult = struct('eventList',eventList,'S_glu',[],'T_glu',[],'dFF_glu',[]);
    return
end

%% ========================= STEP 3: CLUSTER EVENTS → SYNAPSES =========================
% Standard r_min graph connectivity + k-means for large clusters.
% Centroids here already include NMF-derived sub-event centroids from Step 2.
cent = reshape([eventList.centroid_xy], 2, [])';  % Nevent x 2
Nevent = size(cent,1);

% --- Optional UI: tune r_max ---
if opts.doPlot && opts.uiClusterSize && Nevent > 0
    fig = uifigure('Name','Adjust cluster size (r\_max)','Position',[100 100 900 560]);
    ax  = uiaxes(fig,'Position',[20 110 860 420]);
    axis(ax,'ij'); axis(ax,'equal'); xlim(ax,[1 W]); ylim(ax,[1 H]);
    title(ax,'Event centroids colored by synapse cluster (updates with r\_max)');
    xlabel(ax,'x'); ylabel(ax,'y'); hold(ax,'on');
    b = opts.bound;
    rectangle(ax,'Position',[b+1,b+1,W-2*b,H-2*b],'EdgeColor','w','LineWidth',2,'LineStyle','--');
    sld = uislider(fig,'Position',[20 75 860 3],'Limits',[opts.r_max_min opts.r_max_max],'Value',opts.r_max);
    lbl = uilabel(fig,'Position',[20 48 500 22],'Text',sprintf('r_max = %.2f', sld.Value));
    btn = uibutton(fig,'Text','Accept r_max','Position',[700 12 180 32],...
        'ButtonPushedFcn', @(~,~) uiresume(fig));
    fig.UserData.cent  = cent;
    fig.UserData.r_min = opts.r_min;
    fig.UserData.ax    = ax;
    fig.UserData.scatterHandle = [];
    % Wire callback and initial draw AFTER button exists
    sld.ValueChangedFcn = @(src,~) redraw_cluster_plot(fig, lbl, src.Value);
    redraw_cluster_plot(fig, lbl, sld.Value);
    uiwait(fig);
    opts.r_max = sld.Value;
    close(fig);
    fprintf('r_max updated to %.2f\n', opts.r_max);
end

% r_min graph connectivity
r2 = opts.r_min^2;
A  = false(Nevent, Nevent);
for i = 1:Nevent
    dx = cent(i,1) - cent(:,1);
    dy = cent(i,2) - cent(:,2);
    A(i,:) = (dx.^2 + dy.^2) <= r2;
end
G        = graph(A);
synID0   = conncomp(G)';
clusters = accumarray(synID0, (1:Nevent)', [], @(v){v});

synID      = zeros(Nevent,1);
synCounter = 0;

for c = 1:numel(clusters)
    idx = clusters{c};
    pts = cent(idx,:);
    dx  = max(pts(:,1)) - min(pts(:,1));
    dy  = max(pts(:,2)) - min(pts(:,2));
    bboxDiag = hypot(dx, dy);

    if bboxDiag <= opts.r_max || numel(idx) <= 2
        synCounter = synCounter + 1;
        synID(idx) = synCounter;
        continue
    end

    K = max(2, ceil(bboxDiag / opts.r_max));
    try
        lab = kmeans(pts, K, 'Replicates', 5, 'MaxIter', 300);
    catch
        synCounter = synCounter + 1;
        synID(idx) = synCounter;
        continue
    end
    for k = 1:K
        subIdx = idx(lab == k);
        if isempty(subIdx), continue; end
        synCounter = synCounter + 1;
        synID(subIdx) = synCounter;
    end
end

if opts.doPlot && opts.plotCentroids
    fh = figure('Name','Event centroids by synapse ID');
    hold on; axis ij; axis equal; xlim([1 W]); ylim([1 H]);
    b = opts.bound;
    rectangle('Position',[b+1,b+1,W-2*b,H-2*b],'EdgeColor','w','LineWidth',2,'LineStyle','--');
    cmap_syn = hsv(max(synID));
    for s = 1:max(synID)
        ei = find(synID==s);
        if isempty(ei), continue; end
        isSplit = [eventList(ei).nmfSplit];
        scatter(cent(ei(~isSplit),1), cent(ei(~isSplit),2), 12, cmap_syn(s,:), 'filled');
        scatter(cent(ei( isSplit),1), cent(ei( isSplit),2), 20, cmap_syn(s,:), 'd', 'filled');
    end
    title('Event centroids by synapse (diamonds = NMF sub-events)');
    xlabel('x'); ylabel('y'); set(gca,'Color','k'); hold off;
    saveDiag(fh, opts, 'centroids_by_synapse');
end

fprintf('Clustered %d events into %d synapses (r_min=%g, r_max=%g)\n', ...
    Nevent, synCounter, opts.r_min, opts.r_max);


%% ========================= STEP 4: BUILD S_glu0 FROM FOOTPRINT MASKS =========================
%
% CHANGE FROM ORIGINAL:
%   Instead of placing a weighted impulse at the centroid pixel and blurring,
%   we accumulate the actual z-score values (pixVal) from each event's footprint
%   into the corresponding synapse slice of S_glu0.  After summing, we divide
%   by the event count to get a mean z-score footprint, then L2-normalise.
%
S_glu0 = zeros(H, W, synCounter, 'double');
synEventCount0 = zeros(synCounter, 1);

for e = 1:Nevent
    s = synID(e);
    if s == 0, continue; end
    synEventCount0(s) = synEventCount0(s) + 1;

    pixIdx = eventList(e).pixIdx;   % linear indices into H x W
    pixVal = eventList(e).pixVal;   % z-score values (non-binary, >= zThrFoot)

    % Optional: weight by peakZ so high-amplitude events contribute more
    if opts.usePeakWeight
        w = max(eventList(e).peakZ, 0);
        if ~isfinite(w) || w == 0, w = 1; end
    else
        w = 1;
    end

    % Accumulate weighted z-score values into the footprint
    % Accumulate weighted z-score values using safe 2-D slice indexing
    pixIdx_d = double(pixIdx(:));
    pixVal_d = double(pixVal(:));
    fp_slice = zeros(H, W);
    fp_slice(pixIdx_d) = w * pixVal_d;
    S_glu0(:,:,s) = S_glu0(:,:,s) + fp_slice;
        % Equivalent (slower): S_glu0(sub2ind([H W synCounter], r, c, s)) += ...
    end
end

% Divide by event count -> mean weighted z-score footprint per synapse
for s = 1:synCounter
    if synEventCount0(s) > 0
        S_glu0(:,:,s) = S_glu0(:,:,s) / synEventCount0(s);
    end
end

% Optional light smoothing (soft default because footprints are already spatial)
if isfield(opts,'footSigma') && ~isempty(opts.footSigma) && opts.footSigma > 0
    for s = 1:synCounter
        S_glu0(:,:,s) = imgaussfilt(S_glu0(:,:,s), opts.footSigma);
    end
end

% L2 normalise each footprint
for s = 1:synCounter
    v = S_glu0(:,:,s);
    n = sqrt(sum(v(:).^2)) + eps;
    S_glu0(:,:,s) = v / n;
end

% Footprint features for filtering
synArea0 = zeros(synCounter,1);
synAmp0  = zeros(synCounter,1);
for s = 1:synCounter
    fp  = S_glu0(:,:,s);
    thr = max(opts.areaThr, 0);
    synArea0(s) = nnz(fp > thr);
    v = fp(fp > thr);
    synAmp0(s) = ifelse(isempty(v), 0, mean(v));
end

%% ========================= STEP 5: BOUNDARY REMOVAL + FILTERING =========================
bound = opts.bound;
centralSum = squeeze(sum(S_glu0(bound+1:H-bound, bound+1:W-bound, :), [1 2]));
keep = centralSum > 0;
keptIdx = find(keep);

S_glu = S_glu0(:,:,keep);
newNsyn = size(S_glu,3);

mapOldToNew = zeros(synCounter,1);
mapOldToNew(keptIdx) = 1:newNsyn;
newSynID = mapOldToNew(synID);

% Recompute counts / centroids for kept synapses
synEventCount = zeros(newNsyn,1);
synCentroid   = nan(newNsyn,2);
for s = 1:newNsyn
    evIdx = find(newSynID == s);
    synEventCount(s) = numel(evIdx);
    if ~isempty(evIdx)
        synCentroid(s,:) = mean(cent(evIdx,:), 1);
    end
end

% Area / amp for filtering UI
synArea = zeros(newNsyn,1);
synAmp  = zeros(newNsyn,1);
for s = 1:newNsyn
    fp = S_glu(:,:,s);
    nz = fp(fp>0);
    synArea(s) = numel(nz);
    synAmp(s)  = ifelse(isempty(nz), 0, mean(nz));
end

glutrace_tmp = tovec(zMov)' * tovec(S_glu);
amp  = max(glutrace_tmp, [], 1)';
cnt  = synEventCount;
area = synArea;

% Apply saved filter (no UI)
keepUI = true(newNsyn,1);
if ~isempty(opts.synFilter) && ~(opts.doPlot && isfield(opts,'uiSynFilter') && opts.uiSynFilter)
    fr = opts.synFilter;
    keepUI = (amp>=fr.ampRange(1) & amp<=fr.ampRange(2)) & ...
             (cnt>=fr.countRange(1) & cnt<=fr.countRange(2)) & ...
             (area>=fr.areaRange(1) & area<=fr.areaRange(2));
end

% Interactive synapse filter UI
if opts.doPlot && isfield(opts,'uiSynFilter') && opts.uiSynFilter && newNsyn > 0
    fig = uifigure('Name','Synapse filter','Position',[100 100 900 560]);
    ax  = uiaxes(fig,'Position',[20 220 860 310]);
    axis(ax,'ij'); axis(ax,'equal'); xlim(ax,[1 W]); ylim(ax,[1 H]);
    title(ax,'Synapse centroids (red=kept, green=filtered)');
    xlabel(ax,'x'); ylabel(ax,'y'); hold(ax,'on');
    b = opts.bound;
    rectangle(ax,'Position',[b+1,b+1,W-2*b,H-2*b],'EdgeColor','w','LineWidth',2,'LineStyle','--');
    imagesc(ax, zMaxProj > opts.zSeedMin); colormap(ax,gray(256));
    hKeep = scatter(ax, synCentroid(:,1), synCentroid(:,2), 24,'o','MarkerEdgeColor',[1 0 0],'MarkerFaceColor','none');
    hDrop = scatter(ax, synCentroid(:,1), synCentroid(:,2), 24,'o','MarkerEdgeColor',[0 1 0],'MarkerFaceColor','none');
    hold(ax,'off');

    ampL  = fracLimits(amp,  0, 1, eps);
    cntL  = fracLimits(cnt,  0, 1, 1);
    areaL = fracLimits(area, 0, 1, 1);

    [sAmpMin,sAmpMax]   = makePair(fig, 175, 'Max amplitude (fp mean)', ampL, ampL);
    [sCntMin,sCntMax]   = makePair(fig, 130, 'Event count',             cntL, cntL);
    [sAreaMin,sAreaMax] = makePair(fig,  85, 'Area (#fp>0)',           areaL, areaL);
    lbl = uilabel(fig,'Position',[20 50 600 25],'Text','Kept synapses: ...');
    btn = uibutton(fig,'Text','Apply filter','Position',[700 12 180 32],...
        'ButtonPushedFcn', @(~,~) uiresume(fig));

    % Wire callbacks and initial draw AFTER button exists
    sAmpMin.ValueChangedFcn  = @(~,~) redraw2();
    sAmpMax.ValueChangedFcn  = @(~,~) redraw2();
    sCntMin.ValueChangedFcn  = @(~,~) redraw2();
    sCntMax.ValueChangedFcn  = @(~,~) redraw2();
    sAreaMin.ValueChangedFcn = @(~,~) redraw2();
    sAreaMax.ValueChangedFcn = @(~,~) redraw2();
    redraw2();

    uiwait(fig);

    a1=sAmpMin.Value; a2=sAmpMax.Value;
    c1=sCntMin.Value; c2=sCntMax.Value;
    r1=sAreaMin.Value; r2_=sAreaMax.Value;

    opts.synFilter = struct('ampRange',sort([a1 a2]),...
                            'countRange',sort([c1 c2]),...
                            'areaRange',sort([r1 r2_]));
    if opts.saveSynFilter && ~isempty(opts.synFilterFile)
        synFilter = opts.synFilter; %#ok<NASGU>
        save(opts.synFilterFile,'synFilter');
        fprintf('Saved synapse filter to %s\n', opts.synFilterFile);
    end

    keepUI = (amp>=min(a1,a2) & amp<=max(a1,a2)) & ...
             (cnt>=min(c1,c2) & cnt<=max(c1,c2)) & ...
             (area>=min(r1,r2_) & area<=max(r1,r2_));
    close(fig);

    S_glu = S_glu(:,:,keepUI);
    synCentroid   = synCentroid(keepUI,:);
    synEventCount = synEventCount(keepUI);
    newNsyn = size(S_glu,3);

    mapUI = zeros(numel(keepUI),1);
    mapUI(find(keepUI)) = 1:newNsyn;
    valid = newSynID > 0;
    tmp = zeros(size(newSynID));
    tmp(valid) = mapUI(newSynID(valid));
    newSynID = tmp;
end

% Remove synapses with too few events
keepCount = synEventCount >= opts.minEventsSyn;
fprintf('Removing %d synapses with < %d events\n', sum(~keepCount), opts.minEventsSyn);
S_glu = S_glu(:,:,keepCount);
synCentroid   = synCentroid(keepCount,:);
synEventCount = synEventCount(keepCount);
newNsyn = size(S_glu,3);

map2 = zeros(numel(keepCount),1);
map2(find(keepCount)) = 1:newNsyn;
valid = newSynID > 0;
tmp = zeros(size(newSynID));
tmp(valid) = map2(newSynID(valid));
newSynID = tmp;

duration_s        = T / frameRate;
synEventRate_Hz   = synEventCount / max(duration_s, eps);
fprintf('Kept %d synapses. Median event rate: %.3f Hz\n', newNsyn, median(synEventRate_Hz));

%% ========================= STEP 6: VISUALIZE OVERVIEW =========================
RGB = [];
Smap_counts = zeros(H,W,'single');
for e = 1:Nevent
    snew = newSynID(e);
    if snew == 0, continue; end
    Smap_counts(eventList(e).pixIdx) = Smap_counts(eventList(e).pixIdx) + 1;
end

if opts.doPlot
    fh = figure('Name','Event spatial map (counts)');
    imagesc(Smap_counts); axis image tight equal off; colormap hot; colorbar;
    title('Event count map — kept synapses only');
    saveDiag(fh, opts, 'event_count_map');

    cmap_rgb = hsv(newNsyn);
    RGB = zeros(H,W,3,'double');
    for s = 1:newNsyn
        fp = double(S_glu(:,:,s));
        if max(fp(:))==0, continue; end
        V = fp / (max(fp(:)) + eps);
        RGB = max(RGB, bsxfun(@times, reshape(cmap_rgb(s,:),1,1,3), repmat(V,[1 1 3])));
    end

    fh = figure('Name','Synapse footprints (HSV)');
    imshow(RGB); title('Synapse footprints colored by synapse ID (HSV)');
    axis on; hold on;
    validC = ~isnan(synCentroid(:,1));
    plot(synCentroid(validC,1), synCentroid(validC,2),'wo','MarkerSize',7,'LineWidth',0.5);
    hold off;
    saveDiag(fh, opts, 'synapse_footprints_HSV');
end

%% ========================= STEP 7: COMPUTE T_glu =========================
Pixel2compute = find(tovec(max(S_glu,[],3)) > 0);
S_gluVec = tovec(S_glu);
M_gluVec = tovec(double(mov2_sub));
Smat = S_gluVec(Pixel2compute,:);
Msub = M_gluVec(Pixel2compute,:);
Msub(~isfinite(Msub)) = 0;
T_glu = (Smat \ Msub);

F0_glu = []; dFF_glu = [];
if ~isempty(mov2_mc_filt)
    M0_vec = tovec(double(mov2_mc_filt));
    M0_vec(~isfinite(M0_vec)) = 0;
    M0    = M0_vec(Pixel2compute,:);
    F0_glu  = (Smat \ M0);
    dFF_glu = T_glu ./ mean(F0_glu, 2);
end

%% ========================= STEP 8: DIAGNOSTIC PLOTS =========================
% These plots verify correctness:
%   (A) Per-synapse: zoomed footprint side-by-side with its temporal trace
%   (B) Trace comparison: direct pixel projection vs T_glu from LS result
%   (C) Full summary grid: all synapse footprints

if opts.doPlot && newNsyn > 0
    tvec = (0:T-1) / frameRate;  % time axis in seconds

    % ---- (A) Per-synapse zoomed footprint + trace -------------------------
    nShow = min(newNsyn, opts.diagMaxSyn);
    nCols = 4;          % 2 columns: footprint | trace -> 4 actual subplot cols
    nRows = ceil(nShow / 2);
    pad   = opts.diagZoom_pad;

    fh = figure('Name','Per-synapse footprint + trace','Position',[50 50 1400 nRows*220]);
    for s = 1:nShow
        fp   = S_glu(:,:,s);
        xc   = round(synCentroid(s,1));
        yc   = round(synCentroid(s,2));
        xc   = min(max(xc,1),W);
        yc   = min(max(yc,1),H);
        xrng = max(1,xc-pad):min(W,xc+pad);
        yrng = max(1,yc-pad):min(H,yc+pad);
        fpZoom = fp(yrng, xrng);

        row = ceil(s/2);
        col = mod(s-1,2)*2 + 1;  % 1 or 3

        % Footprint subplot
        ax1 = subplot(nRows, nCols, (row-1)*nCols + col);
        imagesc(ax1, fpZoom); axis(ax1,'image'); axis(ax1,'off');
        colormap(ax1, hot); title(ax1, sprintf('Syn %d fp', s), 'FontSize',8);
        % Mark centroid
        hold(ax1,'on');
        plot(ax1, xc - xrng(1)+1, yc - yrng(1)+1, 'c+', 'MarkerSize',8, 'LineWidth',1.5);
        hold(ax1,'off');

        % Trace subplot
        ax2 = subplot(nRows, nCols, (row-1)*nCols + col+1);
        plot(ax2, tvec, T_glu(s,:), 'b', 'LineWidth',0.8);
        xlabel(ax2,'Time (s)','FontSize',7); ylabel(ax2,'AU','FontSize',7);
        title(ax2, sprintf('Syn %d trace (n=%d)', s, synEventCount(s)),'FontSize',8);
        set(ax2,'FontSize',7);
    end
    sgtitle('Per-synapse zoomed footprint (left) & temporal trace (right)','FontSize',11);
    saveDiag(fh, opts, 'per_synapse_footprint_trace');

    % ---- (B) Trace comparison: direct pixel projection vs LS T_glu -------
    % Direct projection: sum(footprint .* movie) / sum(footprint) for each synapse
    % This is the "naive" projection; LS corrects for footprint overlap.
    T_direct = zeros(newNsyn, T);
    for s = 1:newNsyn
        fp   = S_glu(:,:,s);
        fvec = fp(:);
        fsum = sum(fvec.^2) + eps;
        movMat = reshape(double(mov2_sub), H*W, T);
        T_direct(s,:) = (fvec' * movMat) / fsum;
    end

    nComp = min(newNsyn, 6);  % show up to 6 synapses in comparison
    fh = figure('Name','Trace comparison: Direct vs LS','Position',[50 50 1200 nComp*160]);
    for s = 1:nComp
        ax = subplot(nComp, 1, s);
        plot(ax, tvec, T_direct(s,:), 'Color',[0.6 0.6 0.6], 'LineWidth',0.8); hold(ax,'on');
        plot(ax, tvec, T_glu(s,:),    'b', 'LineWidth',1.0); hold(ax,'off');
        ylabel(ax, sprintf('Syn %d',s), 'FontSize',8);
        if s==1
            legend(ax,'Direct projection','LS (T\_glu)','Location','NorthEast','FontSize',7);
            title(ax,'Trace comparison: Direct projection (gray) vs Least-squares T\_glu (blue)');
        end
        set(ax,'FontSize',7);
    end
    xlabel(ax,'Time (s)','FontSize',9);
    saveDiag(fh, opts, 'trace_comparison_direct_vs_LS');

    % ---- (C) Summary grid: all footprints in one figure ------------------
    nc2 = ceil(sqrt(newNsyn));
    nr2 = ceil(newNsyn/nc2);
    fh  = figure('Name','All synapse footprints','Position',[50 50 nc2*120 nr2*120]);
    for s = 1:newNsyn
        subplot(nr2, nc2, s);
        fp   = S_glu(:,:,s);
        xc   = round(synCentroid(s,1));
        yc   = round(synCentroid(s,2));
        xrng = max(1,xc-pad):min(W,xc+pad);
        yrng = max(1,yc-pad):min(H,yc+pad);
        imagesc(fp(yrng,xrng)); axis image off;
        colormap hot;
        title(sprintf('#%d n=%d',s,synEventCount(s)),'FontSize',6);
    end
    sgtitle('All synapse footprints (zoomed, hot colormap)','FontSize',10);
    saveDiag(fh, opts, 'all_synapse_footprints_grid');

    % ---- (D) Event footprints for a random synapse (up to 9 events) ------
    if newNsyn > 0
        sSample = 1;  % show synapse 1 by default (most events)
        evOfSyn = find(newSynID == sSample);
        evOfSyn = evOfSyn(1:min(9,numel(evOfSyn)));
        nc3 = 3; nr3 = ceil(numel(evOfSyn)/3);
        fh = figure('Name', sprintf('Event footprints for Syn %d',sSample),...
                    'Position',[50 50 600 nr3*160]);
        for ii = 1:numel(evOfSyn)
            ev   = eventList(evOfSyn(ii));
            xc   = round(ev.centroid_xy(1));
            yc   = round(ev.centroid_xy(2));
            xrng = max(1,xc-pad):min(W,xc+pad);
            yrng = max(1,yc-pad):min(H,yc+pad);

            zFrame_e = zMov(:,:,ev.frame);
            zZoom    = zFrame_e(yrng,xrng);

            ax = subplot(nr3, nc3, ii);
            imagesc(ax, zZoom); axis(ax,'image'); axis(ax,'off');
            colormap(ax, turbo);
            title(ax, sprintf('Ev %d  t=%d  pZ=%.1f', evOfSyn(ii), ev.frame, ev.peakZ),'FontSize',7);
        end
        sgtitle(sprintf('Individual event footprints (z-score) — Synapse %d', sSample),'FontSize',9);
        saveDiag(fh, opts, sprintf('event_footprints_syn%02d', sSample));
    end

    fprintf('Diagnostic figures saved to: %s\n', opts.diagSaveDir);
end

%% ========================= OUTPUT STRUCT =========================
GluResult = struct();
GluResult.eventList       = eventList;
GluResult.Smap_event      = Smap_counts;
GluResult.S_glu           = S_glu;
GluResult.keep            = keep;
GluResult.synCentroid     = synCentroid;
GluResult.synEventCount   = synEventCount;
GluResult.synEventRate_Hz = synEventRate_Hz;
GluResult.T_glu           = T_glu;
GluResult.F0_glu          = F0_glu;
GluResult.dFF_glu         = dFF_glu;
GluResult.RGB             = RGB;
GluResult.opts            = opts;
GluResult.frameRate       = frameRate;

%% ========================= LOCAL HELPERS =========================
    % ---- Save a figure to diagSaveDir as PNG ----
    function saveDiag(fh, o, tag)
        fname = fullfile(o.diagSaveDir, sprintf('%s_%s.png', o.diagPrefix, tag));
        exportgraphics(fh, fname, 'Resolution', 150);
        fprintf('  Saved: %s\n', fname);
    end

    % ---- ifelse helper ----
    function v = ifelse(cond, a, b)
        if cond; v = a; else; v = b; end
    end

    % ---- Seed-threshold slider redraw ----
    function redraw_slider(sliderHandle)
        zthr = sliderHandle.Value;
        lbl.Text = sprintf('zSeedMin = %.2f', zthr);
        bwMax = max(zMov > zthr, [], 3);
        imagesc(ax, bwMax); colormap(ax,gray);
        axis(ax,'image'); axis(ax,'off'); axis(ax,'tight'); axis(ax,'equal');
        title(ax, sprintf('max(zMov > %.2f, [], 3)', zthr));
        imagesc(ax2, zMaxProj); colormap(ax2,turbo);
        axis(ax2,'image'); axis(ax2,'off'); axis(ax2,'tight'); axis(ax2,'equal');
        title(ax2, 'max(zMov, [], 3)');
        drawnow;
    end

    % ---- Synapse filter UI redraw ----
    function redraw2()
        a1=sAmpMin.Value; a2=sAmpMax.Value;
        c1=sCntMin.Value; c2=sCntMax.Value;
        r1=sAreaMin.Value; r2_=sAreaMax.Value;
        keepUI = (amp>=min(a1,a2) & amp<=max(a1,a2)) & ...
                 (cnt>=min(c1,c2) & cnt<=max(c1,c2)) & ...
                 (area>=min(r1,r2_) & area<=max(r1,r2_));
        x = synCentroid(:,1); y = synCentroid(:,2);
        hKeep.XData = x(keepUI);  hKeep.YData = y(keepUI);
        hDrop.XData = x(~keepUI); hDrop.YData = y(~keepUI);
        lbl.Text = sprintf('Kept synapses: %d / %d', sum(keepUI), numel(keepUI));
        drawnow limitrate;
    end

    % ---- Slider pair helper ----
    function [sMin,sMax] = makePair(figH, y, label, lims, v0)
        uilabel(figH,'Text',label,'Position',[20 y 160 20]);
        sMin = uislider(figH,'Position',[180 y+8 300 3],'Limits',lims,'Value',v0(1));
        sMax = uislider(figH,'Position',[520 y+8 300 3],'Limits',lims,'Value',v0(2));
        uilabel(figH,'Text','min','Position',[180 y-12 30 20]);
        uilabel(figH,'Text','max','Position',[520 y-12 30 20]);
    end

    % ---- Cluster-size slider redraw ----
    function redraw_cluster_plot(figH, lblH, r_max_val)
        lblH.Text = sprintf('r_max = %.2f', r_max_val);
        centL  = figH.UserData.cent;
        r_minL = figH.UserData.r_min;
        axL    = figH.UserData.ax;
        synID_tmp = cluster_by_rmin_rmax(centL, r_minL, r_max_val);
        K = max(synID_tmp);
        cmap2 = hsv(max(K,1));
        cmap2 = cmap2(randi(size(cmap2 ,1),size(cmap2 ,1),1),:);

        C = cmap2(synID_tmp,:);
        if isempty(figH.UserData.scatterHandle) || ~isgraphics(figH.UserData.scatterHandle)
            figH.UserData.scatterHandle = scatter(axL, centL(:,1), centL(:,2), 12, C,'filled');
        else
            h = figH.UserData.scatterHandle;
            h.XData = centL(:,1); h.YData = centL(:,2); h.CData = C;
        end
        title(axL, sprintf('Event centroids (r_{max}=%.2f, #clusters=%d)', r_max_val, K));
        drawnow;
    end

    % ---- Centroid clustering with r_min + r_max (used by UI only) ----
    function synIDout = cluster_by_rmin_rmax(centIn, r_minIn, r_maxIn)
        N = size(centIn,1);
        r2 = r_minIn^2;
        Ain = false(N,N);
        for i = 1:N
            dxi = centIn(i,1)-centIn(:,1); dyi = centIn(i,2)-centIn(:,2);
            Ain(i,:) = (dxi.^2+dyi.^2) <= r2;
        end
        G2 = graph(Ain);
        sID0 = conncomp(G2)';
        cls  = accumarray(sID0,(1:N)',[], @(v){v});
        synIDout = zeros(N,1);
        sc2 = 0;
        for ci = 1:numel(cls)
            idxI = cls{ci}; ptsI = centIn(idxI,:);
            dxI = max(ptsI(:,1))-min(ptsI(:,1)); dyI = max(ptsI(:,2))-min(ptsI(:,2));
            bd  = hypot(dxI,dyI);
            if bd<=r_maxIn || numel(idxI)<=2
                sc2=sc2+1; synIDout(idxI)=sc2; continue
            end
            KI = max(2,ceil(bd/r_maxIn));
            try
                labI = kmeans(ptsI,KI,'Replicates',3,'MaxIter',200);
            catch
                sc2=sc2+1; synIDout(idxI)=sc2; continue
            end
            for kk=1:KI
                sI2=idxI(labI==kk); if isempty(sI2),continue;end
                sc2=sc2+1; synIDout(sI2)=sc2;
            end
        end
    end

    % ---- Percentile-based slider limits ----
    function L = fracLimits(v, fLo, fHi, epsPad)
        v = v(isfinite(v)); v = sort(v(:)); n = numel(v);
        if n==0, L=[0 1]; return; end
        iLo = max(1,min(n,floor(fLo*n)));
        iHi = max(1,min(n,ceil(fHi*n)));
        if iHi<iLo, iHi=iLo; end
        L = [v(iLo) v(iHi)];
        if L(1)==L(2), L=[L(1) L(1)+epsPad]; end
    end
end  % extractGluSNFR2