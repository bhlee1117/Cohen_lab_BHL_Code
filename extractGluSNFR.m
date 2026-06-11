function GluResult = extractGluSNFR(mov2_sub, mov2_mc_filt, exposuretime2, opts, S_glu_in)
% extractGluSNFR
%   Detect pixel-wise iGluSNFR events, cluster into synapses, build S_glu (H×W×Nsyn),
%   remove boundary synapses, visualize, and compute T_glu by least-squares projection.
%
% USAGE
%   GluResult = extractGluSNFR(mov2_sub, mov2_mc_filt, exposuretime2);
%   GluResult = extractGluSNFR(mov2_sub, mov2_mc_filt, exposuretime2, opts);
%
% INPUTS
%   mov2_sub      : H×W×T movie (e.g., background-subtracted / residual movie)
%   mov2_mc_filt  : H×W×T movie used for baseline (F0) projection (can be [])
%   exposuretime2 : seconds per frame
%   opts          : struct of parameters (optional)
%
% OUTPUT (struct GluResult)
%   .eventList        : detected event regions (frame, pixIdx, centroid_xy, peakZ)
%   .Smap_event       : H×W count map of events (kept synapses only)
%   .S_glu            : H×W×Nsyn normalized synapse footprints (after boundary removal)
%   .synCentroid      : Nsyn×2 centroid positions [x y]
%   .synEventCount    : Nsyn×1 number of events per synapse
%   .synEventRate_Hz  : Nsyn×1 event rates  
%   .keep             : keep mask (before remap)
%   .T_glu            : Nsyn×T temporal projection traces
%   .F0_glu           : Nsyn×T baseline projection traces (empty if mov2_mc_filt==[])
%   .dFF_glu          : Nsyn×T dF/F (empty if mov2_mc_filt==[])
%   .RGB              : H×W×3 HSV-colored footprint visualization (kept only)

%% ---------------- defaults ----------------
if nargin < 4 || isempty(opts), opts = struct; end
    function val = getOpt(opts, field, default)
        if isfield(opts, field) && ~isempty(opts.(field))
            val = opts.(field);
        else
            val = default;
        end
    end

frameRate     = 1/exposuretime2;

% z-score / peak seed params
opts.smoothSigma   = getOpt(opts,'smoothSigma',   1);
opts.zSeedMin      = getOpt(opts,'zSeedMin',      3);
opts.zThrFoot      = getOpt(opts,'zThrFoot',      2.5);
opts.minPeakDist_s = getOpt(opts,'minPeakDist_s', 0.2);
opts.minPeakProm   = getOpt(opts,'minPeakProm',   0.75);
opts.minArea       = getOpt(opts,'minArea',       2);
opts.eventPixMax   = getOpt(opts,'eventPixMax',   5);
opts.r_min         = getOpt(opts,'r_min',         3);
opts.r_max         = getOpt(opts,'r_max',         8);
opts.bound         = getOpt(opts,'bound',         10);
opts.doPlot         = getOpt(opts,'doPlot', true);   % existing
opts.uiSeedThresh   = getOpt(opts,'uiSeedThresh', true);  % (1) slider for zSeedMin
opts.uiSynFilter    = getOpt(opts,'uiSynFilter', true);   % (4) interactive synapse filter
opts.plotCentroids  = getOpt(opts,'plotCentroids', true); % (3) centroid scatter
opts.keepPixWeights = getOpt(opts,'keepPixWeights', true); % (2) keep pixel values as weights
opts.usePeakWeight = getOpt(opts,'usePeakWeight', true);
opts.minEventsSyn  = getOpt(opts,'minEventsSyn', 3);   % keep synapses with >= this many events

opts.uiClusterSize = getOpt(opts,'uiClusterSize', true);  % new
opts.r_max_min     = getOpt(opts,'r_max_min',  0.001);
opts.r_max_max     = getOpt(opts,'r_max_max', 10);

opts.footSigma = getOpt(opts,'footSigma', 1);  % footprint blur of centroid impulses
opts.areaThr   = getOpt(opts,'areaThr',   0);    % if blurred, set e.g. 1e-6 or prctile-based
opts.synFilter = getOpt(opts,'synFilter', []);   % stored thresholds
opts.saveSynFilter = getOpt(opts,'saveSynFilter', false);
opts.synFilterFile = getOpt(opts,'synFilterFile', '');


%%
if ~isfield(opts,'mode') || isempty(opts.mode)
    opts.mode = 'detect';  % 'detect' or 'project'
end

if strcmpi(opts.mode, 'project')
    if nargin < 5 || isempty(S_glu_in)
        error('Project mode requires S_glu_in (H x W x Nsyn).');
    end

    % --- Project-only: compute T_glu, F0_glu, dFF_glu ---
    mov = double(mov2_sub);
    [H,W,T] = size(mov);

    % ---- Validate movie dimensions ----
    if ~isempty(mov2_mc_filt)
        sz1 = size(mov2_sub);
        sz2 = size(mov2_mc_filt);

        if numel(sz1) ~= 3 || numel(sz2) ~= 3
            error('mov2_sub and mov2_mc_filt must both be H×W×T 3D arrays.');
        end

        if any(sz1 ~= sz2)
            error(['Size mismatch between mov2_sub and mov2_mc_filt.\n' ...
                'mov2_sub size:      [%d %d %d]\n' ...
                'mov2_mc_filt size:  [%d %d %d]'], ...
                sz1(1),sz1(2),sz1(3), ...
                sz2(1),sz2(2),sz2(3));
        end
    end


    [Hs,Ws,Nsyn] = size(S_glu_in);
    assert(H==Hs && W==Ws, 'S_glu_in must match movie spatial size.');

    Npix = H*W;
    S = reshape(double(S_glu_in), Npix, Nsyn);  % Npix x Nsyn
    M = reshape(mov, Npix, T);                  % Npix x T
    M(~isfinite(M)) = 0;

    % Least squares projection (stable QR)
    T_glu = (S \ M);                            % Nsyn x T

    F0_glu = [];
    dFF_glu = [];
    if ~isempty(mov2_mc_filt)
        M0 = reshape(double(mov2_mc_filt), Npix, T);
        M0(~isfinite(M0)) = 0;
        F0_glu = (S \ M0);
        dFF_glu = T_glu ./ mean(F0_glu, 2);
    end

    GluResult = struct();
    GluResult.S_glu = S_glu_in;
    GluResult.T_glu = T_glu;
    GluResult.F0_glu = F0_glu;
    GluResult.dFF_glu = dFF_glu;
    GluResult.frameRate = 1/exposuretime2;
    GluResult.opts = opts;
    return
end
%% ---------------- prep ----------------
mov = double(mov2_sub);

[H,W,T] = size(mov);
minPeakDist = max(1, round(opts.minPeakDist_s * frameRate));

mov = imgaussfilt(mov,opts.smoothSigma);

% Robust per-pixel z-score movie
medImg = median(mov, 3);
madImg = mad(mov, 1, 3);
sigmaImg = 1.4826 * madImg + eps;
%sigmaImg = imgaussfilt(mean(mov2_mc_filt,3),2);
zMov = (mov - medImg) ./ sigmaImg;
% Precompute context projection
    zMaxProj = max(zMov, [], 3, 'omitnan');
    zMaxProj(zMaxProj==0)=median(tovec(zMaxProj(zMaxProj~=0)));

if opts.doPlot && opts.uiSeedThresh

    zSeedMin0 = opts.zSeedMin;
    zMin = 0;
    zMax = max(6, ceil(max(zMov(:),[],'omitnan')));

    fig = uifigure('Name','Adjust zSeedMin','Position',[100 100 900 450]);
    ax  = uiaxes(fig,'Position',[20 60 420 360]);
    ax2 = uiaxes(fig,'Position',[460 60 420 360]);

    % Create slider
    sld = uislider(fig,...
        'Position',[20 30 860 3],...
        'Limits',[zMin zMax],...
        'Value',zSeedMin0);

    lbl = uilabel(fig,...
        'Position',[20 5 400 20],...
        'Text',sprintf('zSeedMin = %.2f', sld.Value));

    % Attach callback safely
    sld.ValueChangedFcn = @(src,event) redraw_slider(src);

    % Initial draw
    redraw_slider(sld);

    % Accept button
    btn = uibutton(fig,...
        'Text','Use this zSeedMin',...
        'Position',[740 5 140 25],...
        'ButtonPushedFcn', @(btn,event) uiresume(fig));

    uiwait(fig);

    opts.zSeedMin = sld.Value;
    close(fig);

    fprintf('zSeedMin updated to %.2f\n', opts.zSeedMin);
end

fprintf('Finding per-pixel peaks (this can take a bit)...\n');

%% ---------------- 1) collect peak seeds by frame ----------------
seedByFrame = cell(T,1);

for y = 1:H
    for x = 1:W
        tr = squeeze(zMov(y,x,:));
        if max(tr) < opts.zSeedMin
            continue
        end

        [pks, locs] = findpeaks(tr, ...
            'MinPeakDistance',   minPeakDist, ...
            'MinPeakProminence', opts.minPeakProm);

        if isempty(locs), continue; end
        keepPk = pks >= opts.zSeedMin;
        locs = locs(keepPk);
        if isempty(locs), continue; end

        lin = sub2ind([H,W], y, x);
        for k = 1:numel(locs)
            seedByFrame{locs(k)}(end+1) = lin; %#ok<AGROW>
        end
    end
end

nSeedsTotal = sum(cellfun(@numel, seedByFrame));
fprintf('Total seed peaks found: %d\n', nSeedsTotal);

%% ---------------- 2) segment footprints at seed frames ----------------
eventList = struct('frame', {}, 'pixIdx', {}, 'centroid_xy', {}, 'peakZ', {});
compCount = 0;

for t = 1:T
    seeds = seedByFrame{t};
    if isempty(seeds), continue; end

    zFrame = zMov(:,:,t);

    bw = zFrame > opts.zThrFoot;
    bw = bwareaopen(bw, opts.minArea);
    bw = imopen(bw, strel('disk', 1));
    bw = imfill(bw, 'holes');

    CC = bwconncomp(bw, 8);
    if CC.NumObjects == 0, continue; end

    L = labelmatrix(CC);  % HxW
    seedLabels = unique(L(seeds));
    seedLabels(seedLabels==0) = [];
    if isempty(seedLabels), continue; end

    for lab = seedLabels(:)'
        pixIdx = find(L == lab);
        nPix = numel(pixIdx);
        if nPix < opts.minArea, continue; end
        [r,c] = ind2sub([H,W], pixIdx);
        compCount = compCount + 1;
        eventList(compCount).frame = t;
        eventList(compCount).pixIdx = pixIdx;
        eventList(compCount).centroid_xy = [mean(c) mean(r)];
        eventList(compCount).peakZ = max(zFrame(pixIdx));
        eventList(compCount).pixVal = zFrame(pixIdx);   % if you're keeping pixel weights
        eventList(compCount).meanZ  = mean(zFrame(pixIdx));
    end
end

fprintf('Extracted %d event regions from pixel-trace peaks\n', numel(eventList));

if isempty(eventList)
    warning('No events detected. Returning empty results.');
    GluResult = struct('eventList',eventList,'S_glu',[],'T_glu',[],'dFF_glu',[]);
    return
end

%% ---------------- 3) cluster event centroids into synapses (r_min + r_max) ----------------
cent = reshape([eventList.centroid_xy], 2, [])';   % Nevent x 2
Nevent = size(cent,1);

if opts.doPlot && opts.uiClusterSize && Nevent > 0
    % --- interactive r_max tuning (max synapse "width") ---
    fig = uifigure('Name','Adjust cluster size (r\_max)','Position',[100 100 900 520]);
    ax  = uiaxes(fig,'Position',[20 80 860 420]);
    axis(ax,'ij'); axis(ax,'equal');
    xlim(ax,[1 W]); ylim(ax,[1 H]);
    title(ax,'Event centroids colored by synapse cluster (updates with r\_max)');
    xlabel(ax,'x'); ylabel(ax,'y'); hold(ax,'on');

    % boundary box
    b = opts.bound;
    rectangle(ax,'Position',[b+1, b+1, W-2*b, H-2*b], ...
        'EdgeColor','w','LineWidth',2,'LineStyle','--');

    % slider + label
    sld = uislider(fig, ...
        'Position',[20 50 860 3], ...
        'Limits',[opts.r_max_min opts.r_max_max], ...
        'Value', opts.r_max);
    lbl = uilabel(fig,'Position',[20 20 500 22], ...
        'Text', sprintf('r_max = %.2f', sld.Value));

    % store data in fig.UserData (avoid scoping problems)
    fig.UserData.cent = cent;
    fig.UserData.r_min = opts.r_min;
    fig.UserData.ax = ax;
    fig.UserData.scatterHandle = [];  % will hold one scatter (fast redraw)

    % initial draw
    redraw_cluster_plot(fig, lbl, sld.Value);

    % update on slider
    sld.ValueChangedFcn = @(src,evt) redraw_cluster_plot(fig, lbl, src.Value);

    % accept button
    uibutton(fig,'Text','Use this r_max','Position',[760 15 120 30], ...
        'ButtonPushedFcn', @(btn,event) uiresume(fig));

    uiwait(fig);

    opts.r_max = sld.Value;     % <-- this is the chosen value used downstream
    close(fig);
    fprintf('r_max updated to %.2f\n', opts.r_max);
end


r2 = opts.r_min^2;
A = false(Nevent, Nevent);
for i = 1:Nevent
    dx = cent(i,1) - cent(:,1);
    dy = cent(i,2) - cent(:,2);
    A(i,:) = (dx.^2 + dy.^2) <= r2;
end
G = graph(A);
synID0 = conncomp(G)';  % initial
clusters = accumarray(synID0, (1:Nevent)', [], @(v){v});

synID = zeros(Nevent,1);
synCounter = 0;

for c = 1:numel(clusters)
    idx = clusters{c};
    pts = cent(idx,:);

    dx = max(pts(:,1)) - min(pts(:,1));
    dy = max(pts(:,2)) - min(pts(:,2));
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
        subPts = cent(subIdx,:);

        dx = max(subPts(:,1)) - min(subPts(:,1));
        dy = max(subPts(:,2)) - min(subPts(:,2));
        subDiag = hypot(dx, dy);

        if subDiag <= opts.r_max || numel(subIdx) <= 2
            synCounter = synCounter + 1;
            synID(subIdx) = synCounter;
        else
            K2 = max(2, ceil(subDiag / opts.r_max));
            lab2 = kmeans(subPts, K2, 'Replicates', 5, 'MaxIter', 300);
            for kk = 1:K2
                subsubIdx = subIdx(lab2 == kk);
                if isempty(subsubIdx), continue; end
                synCounter = synCounter + 1;
                synID(subsubIdx) = synCounter;
            end
        end
    end
end

if opts.doPlot && opts.plotCentroids
    figure; hold on; axis ij; axis equal;
    xlim([1 W]); ylim([1 H]);

    % Draw boundary box for allowed central region
    b = opts.bound;
    rectangle('Position',[b+1, b+1, W-2*b, H-2*b], ...
        'EdgeColor','w','LineWidth',2,'LineStyle','--');

    % Color by synapse cluster
    cmap = hsv(max(synID));
    for s = 1:max(synID)
        idx = find(synID==s);
        if isempty(idx), continue; end
        scatter(cent(idx,1), cent(idx,2), 12, cmap(s,:), 'filled');
    end

    title('Event centroids colored by clustered synapse ID');
    xlabel('x'); ylabel('y');
    set(gca,'Color','k');  % optional: dark background
    hold off;
end

fprintf('Clustered %d events into %d synapses (r_min=%g, r_max=%g)\n', ...
    Nevent, synCounter, opts.r_min, opts.r_max);

%% ---------------- 4) build S_glu (H×W×Nsyn) from event centroids only ----------------
S_glu0 = zeros(H, W, synCounter, 'double');
synEventCount0 = zeros(synCounter,1);

for e = 1:Nevent
    s = synID(e);
    synEventCount0(s) = synEventCount0(s) + 1;

    % --- centroid (xy) -> nearest pixel ---
    xy = eventList(e).centroid_xy;   % [x y]
    xC = round(xy(1));
    yC = round(xy(2));

    % clamp to image bounds
    xC = min(max(xC,1), W);
    yC = min(max(yC,1), H);

    % --- event weight ---
    if opts.usePeakWeight && isfield(eventList,'peakZ')
        w = max(eventList(e).peakZ, 0);
        if ~isfinite(w) || w==0, w = 1; end
    else
        w = 1;
    end

    % add weight only at centroid pixel
    S_glu0(yC, xC, s) = S_glu0(yC, xC, s) + w;
end

% average by event count (so footprint is "mean centroid weight per event")
for s = 1:synCounter
    if synEventCount0(s) > 0
        S_glu0(:,:,s) = S_glu0(:,:,s) / synEventCount0(s);
    end
end

% OPTIONAL: spatial smoothing of centroid impulse map to form a blob footprint
% (recommended; otherwise footprint is very sparse and LS projection is noisy)
if isfield(opts,'footSigma') && ~isempty(opts.footSigma) && opts.footSigma > 0
    for s = 1:synCounter
        S_glu0(:,:,s) = imgaussfilt(S_glu0(:,:,s), opts.footSigma);
    end
end

% L2 normalize each synapse footprint
for s = 1:synCounter
    v = S_glu0(:,:,s);
    n = sqrt(sum(v(:).^2)) + eps;
    S_glu0(:,:,s) = v / n;
end

% ---- features for filtering ----
synArea0 = zeros(synCounter,1);
synAmp0  = zeros(synCounter,1);

for s = 1:synCounter
    fp = S_glu0(:,:,s);

    % If you smoothed, area thresholding by >0 becomes meaningless.
    % Use a small percentile/threshold instead:
    if isfield(opts,'areaThr') && ~isempty(opts.areaThr)
        thr = opts.areaThr;
    else
        thr = 0;  % if no smoothing, impulses are exact zeros elsewhere
    end

    synArea0(s) = nnz(fp > thr);

    v = fp(fp > thr);
    if isempty(v), synAmp0(s) = 0;
    else,          synAmp0(s) = mean(v);
    end
end

%% ---------------- 5) remove boundary synapses + optional interactive synapse filtering ----------------
bound = opts.bound;

% ---- 5a) boundary keep (unchanged logic) ----
centralSum = squeeze(sum(S_glu0(bound+1:H-bound, bound+1:W-bound, :), [1 2]));
keep = centralSum > 0;                 % keep mask over synCounter synapses

keptIdx = find(keep);                  % original synapse indices after boundary filter
S_glu = S_glu0(:,:,keep);
newNsyn = size(S_glu,3);

% remap old synapse index -> new synapse index for events
mapOldToNew = zeros(synCounter,1);
mapOldToNew(keptIdx) = 1:newNsyn;
newSynID = mapOldToNew(synID);         % event->(kept synapse id), 0 means dropped by boundary

% recompute synEventCount, synCentroid for kept synapses
synEventCount = zeros(newNsyn,1);
synCentroid = nan(newNsyn,2);
for s = 1:newNsyn
    evIdx = find(newSynID == s);
    synEventCount(s) = numel(evIdx);
    if ~isempty(evIdx)
        synCentroid(s,:) = mean(cent(evIdx,:), 1);
    end
end

% ---- 5a-2) compute synapse features for interactive filtering (kept only) ----
% area: count of pixels with footprint > 0
synArea = zeros(newNsyn,1);
% amp: mean footprint value over nonzero pixels (proxy for mean amplitude)
synAmp  = zeros(newNsyn,1);

for s = 1:newNsyn
    fp = S_glu(:,:,s);
    nz = fp(fp>0);
    synArea(s) = numel(nz);
    if isempty(nz)
        synAmp(s) = 0;
    else
        synAmp(s) = mean(nz);
    end
end

glutrace_tmp=tovec(zMov)'*tovec(S_glu);
amp  = max(glutrace_tmp,[],1)';
cnt  = synEventCount;
area = synArea;

% ---- 5a-3) OPTIONAL interactive filtering by [amp, count, area] ----
% This produces a keepUI mask over the *kept synapses* (1..newNsyn).
keepUI = true(newNsyn,1);

% ---- If synFilter exists and UI is OFF, apply saved thresholds ----
if ~isempty(opts.synFilter) && ...
        ~(opts.doPlot && isfield(opts,'uiSynFilter') && opts.uiSynFilter)

    fr = opts.synFilter;
    keepUI = (amp>=fr.ampRange(1) & amp<=fr.ampRange(2)) & ...
        (cnt>=fr.countRange(1) & cnt<=fr.countRange(2)) & ...
        (area>=fr.areaRange(1) & area<=fr.areaRange(2));
end

if opts.doPlot && isfield(opts,'uiSynFilter') && opts.uiSynFilter && newNsyn > 0
    fig = uifigure('Name','Synapse filter','Position',[100 100 900 520]);
    % ---- main axis: centroid map colored by synapse ----
    ax = uiaxes(fig,'Position',[20 200 860 300]);
    axis(ax,'ij'); axis(ax,'equal');
    xlim(ax,[1 W]); ylim(ax,[1 H]);
    title(ax,'Synapse centroids (colored by synapse; gray = filtered out)');
    xlabel(ax,'x'); ylabel(ax,'y');
    hold(ax,'on');

    % Draw boundary box (central keep region)
    b = opts.bound;
    rectangle(ax,'Position',[b+1, b+1, W-2*b, H-2*b], ...
        'EdgeColor','w','LineWidth',2,'LineStyle','--');

    % Build color map per synapse
    cmap = hsv(newNsyn);

    % Pre-create one scatter object per synapse (fast updates later)
    imagesc(ax,zMaxProj>opts.zSeedMin); hold all;
    colormap(ax,gray(256))

    % Create two scatter handles ONCE
    hKeep = scatter(ax, synCentroid(:,1), synCentroid(:,2), 24, 'o', ...
        'MarkerEdgeColor',[1 0 0], 'MarkerFaceColor','none');  % kept (red)
    hold(ax,'on');
    hDrop = scatter(ax, synCentroid(:,1), synCentroid(:,2), 24, 'o', ...
        'MarkerEdgeColor',[0 1 0], 'MarkerFaceColor','none');  % filtered (green)
    hold(ax,'off');

    % --- Narrow slider limits using linear fractions of sorted values ---

    ampL  = fracLimits(amp,  0, 1, eps);
    cntL  = fracLimits(cnt,  0, 1, 1);
    areaL = fracLimits(area, 0, 1, 1);

    % --- Initial slider values: start at the narrowed limits ---
    amp0  = ampL;
    cnt0  = cntL;
    area0 = areaL;

    [sAmpMin,sAmpMax]   = makePair(160,'Max amplitude (fp mean)', ampL,  amp0);
    [sCntMin,sCntMax]   = makePair(120,'Event count',             cntL,  cnt0);
    [sAreaMin,sAreaMax] = makePair( 80,'Area (#fp>0)',            areaL, area0);

    lbl = uilabel(fig,'Position',[20 30 600 25],'Text','Kept synapses: ...');

    % Wire callbacks
    sAmpMin.ValueChangedFcn   = @(~,~) redraw2();
    sAmpMax.ValueChangedFcn   = @(~,~) redraw2();
    sCntMin.ValueChangedFcn   = @(~,~) redraw2();
    sCntMax.ValueChangedFcn   = @(~,~) redraw2();
    sAreaMin.ValueChangedFcn  = @(~,~) redraw2();
    sAreaMax.ValueChangedFcn  = @(~,~) redraw2();

    redraw2();

    uibutton(fig,'Text','Apply filter','Position',[760 30 120 30], ...
        'ButtonPushedFcn', @(btn,event) uiresume(fig));

    uiwait(fig);

    % Final keepUI based on slider values at acceptance
    a1=sAmpMin.Value; a2=sAmpMax.Value;
    c1=sCntMin.Value; c2=sCntMax.Value;
    r1=sAreaMin.Value; r2=sAreaMax.Value;

    % Store thresholds in a struct
    opts.synFilter = struct();
    opts.synFilter.ampRange  = sort([a1 a2]);
    opts.synFilter.countRange = sort([c1 c2]);
    opts.synFilter.areaRange  = sort([r1 r2]);

    % Optionally save to disk
    if opts.saveSynFilter && ~isempty(opts.synFilterFile)
        synFilter = opts.synFilter; %#ok<NASGU>
        save(opts.synFilterFile,'synFilter');
        fprintf('Saved synapse filter to %s\n', opts.synFilterFile);
    end


    keepUI = (amp>=min(a1,a2) & amp<=max(a1,a2)) & ...
        (cnt>=min(c1,c2) & cnt<=max(c1,c2)) & ...
        (area>=min(r1,r2) & area<=max(r1,r2));

    close(fig);

    % Apply UI mask to kept synapses: shrink S_glu etc.
    S_glu = S_glu(:,:,keepUI);
    synCentroid = synCentroid(keepUI,:);
    synEventCount = synEventCount(keepUI);

    % Update newNsyn after UI filtering
    newNsyn = size(S_glu,3);

    % Update newSynID (events) to reflect removed synapses:
    % newSynID currently refers to indices 1..(old newNsyn before keepUI).
    mapUI = zeros(numel(keepUI),1);
    mapUI(find(keepUI)) = 1:newNsyn;

    valid = newSynID > 0;
    tmp = zeros(size(newSynID));
    tmp(valid) = mapUI(newSynID(valid));
    newSynID = tmp;
end

% ---- 5b) remove synapses with too few events (your existing credibility filter) ----
keepCount = synEventCount >= opts.minEventsSyn;

fprintf('Removing %d synapses with < %d events\n', sum(~keepCount), opts.minEventsSyn);

S_glu = S_glu(:,:,keepCount);
synCentroid = synCentroid(keepCount,:);
synEventCount = synEventCount(keepCount);

newNsyn = size(S_glu,3);

% Update newSynID mapping again after keepCount
map2 = zeros(numel(keepCount),1);
map2(find(keepCount)) = 1:newNsyn;
valid = newSynID > 0;
tmp = zeros(size(newSynID));
tmp(valid) = map2(newSynID(valid));
newSynID = tmp;

duration_s = T / frameRate;
synEventRate_Hz = synEventCount / max(duration_s, eps);

fprintf('Kept %d synapses (removed %d). Median event rate: %.3f Hz\n', ...
    newNsyn, synCounter-newNsyn, median(synEventRate_Hz));

%% ---------------- 6) visualization (optional) ----------------
RGB = [];
Smap_counts = zeros(H,W,'single');

for e = 1:Nevent
    snew = newSynID(e);
    if snew == 0, continue; end
    Smap_counts(eventList(e).pixIdx) = Smap_counts(eventList(e).pixIdx) + 1;
end

if opts.doPlot
    figure; imagesc(Smap_counts); axis image tight equal off; colormap hot; colorbar;
    title('Event spatial map (counts) — kept synapses only');

    cmap = hsv(newNsyn);
    RGB = zeros(H, W, 3, 'double');

    for s = 1:newNsyn
        fp = double(S_glu(:,:,s));
        if max(fp(:)) == 0, continue; end
        V = fp / (max(fp(:)) + eps);
        color_rgb = reshape(cmap(s,:), 1,1,3);
        RGB_syn = bsxfun(@times, color_rgb, repmat(V,[1,1,3]));
        RGB = max(RGB, RGB_syn);
    end

    figure; imshow(RGB);
    title('Normalized Synapse Footprints colored by hsv(Nsyn) (kept only)');
    axis on; hold on;
    validC = ~isnan(synCentroid(:,1));
    plot(synCentroid(validC,1), synCentroid(validC,2), 'wo', 'MarkerSize', 7, 'LineWidth', 0.5);
    hold off;
end

%% ---------------- 7) compute T_glu by least squares projection ----------------
% NOTE: do NOT use tovec(S_glu)\tovec(mov) (dimension mismatch).
% Correct: reshape to (Npix x Nsyn) and (Npix x T), then solve S\Movie.
Pixel2compute=find(tovec(max(S_glu,[],3))>0);
Npix = length(Pixel2compute);
S_gluVec=tovec(S_glu);
M_gluVec=tovec(double(mov2_sub));
Smat =S_gluVec(Pixel2compute,:);
Msub =M_gluVec(Pixel2compute,:);

Msub(~isfinite(Msub)) = 0;
T_glu = (Smat \ Msub);

F0_glu = [];
dFF_glu = [];
if ~isempty(mov2_mc_filt)
    M0_vec = tovec(double(mov2_mc_filt));
    M0_vec(~isfinite(M0_vec)) = 0;
    M0 = M0_vec(Pixel2compute,:);
    F0_glu = (Smat \ M0);
    dFF_glu = T_glu ./ mean(F0_glu, 2);
end

% Npix = H*W;
% Smat = reshape(S_glu, Npix, newNsyn);      % Npix x Nsyn
% Msub = reshape(double(mov2_sub), Npix, T); % Npix x T
% 
% % Remove NaN/Inf just in case
% 
% 
% % Least squares (MATLAB uses QR). If ill-conditioned, add ridge outside this function.
% T_glu = (Smat \ Msub);                     % Nsyn x T
% 
% F0_glu = [];
% dFF_glu = [];
% if ~isempty(mov2_mc_filt)
%     M0 = reshape(double(mov2_mc_filt), Npix, T);
%     M0(~isfinite(M0)) = 0;
%     F0_glu = (Smat \ M0);
%     dFF_glu = T_glu ./ mean(F0_glu, 2);
% end

%% ---------------- outputs ----------------
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

% ---- SAFE redraw using passed slider handle ----
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

    function redraw2()
    a1=sAmpMin.Value; a2=sAmpMax.Value;
    c1=sCntMin.Value; c2=sCntMax.Value;
    r1=sAreaMin.Value; r2=sAreaMax.Value;

    keepUI = (amp>=min(a1,a2) & amp<=max(a1,a2)) & ...
             (synEventCount>=min(c1,c2) & synEventCount<=max(c1,c2)) & ...
             (synArea>=min(r1,r2) & synArea<=max(r1,r2));

    x = synCentroid(:,1);
    y = synCentroid(:,2);

    % Update kept scatter
    hKeep.XData = x(keepUI);
    hKeep.YData = y(keepUI);

    % Update dropped scatter (either hide or show)
    if isfield(opts,'uiHideFiltered') && opts.uiHideFiltered
        hDrop.Visible = 'off';
    else
        hDrop.Visible = 'on';
        hDrop.XData = x(~keepUI);
        hDrop.YData = y(~keepUI);
    end

    lbl.Text = sprintf('Kept synapses: %d / %d', sum(keepUI), numel(keepUI));
    drawnow limitrate;   % faster than drawnow
end


    function [sMin,sMax] = makePair(y, label, lims, v0)
        uilabel(fig,'Text',label,'Position',[20 y 160 20]);
        sMin = uislider(fig,'Position',[180 y+8 300 3], 'Limits',lims, 'Value',v0(1));
        sMax = uislider(fig,'Position',[520 y+8 300 3], 'Limits',lims, 'Value',v0(2));
        uilabel(fig,'Text','min','Position',[180 y-12 30 20]);
        uilabel(fig,'Text','max','Position',[520 y-12 30 20]);
    end

    function redraw_cluster_plot(fig, lbl, r_max_val)
        lbl.Text = sprintf('r_max = %.2f', r_max_val);

        cent = fig.UserData.cent;
        r_min = fig.UserData.r_min;
        ax = fig.UserData.ax;

        % Recompute clustering with current r_max
        synID_tmp = cluster_by_rmin_rmax(cent, r_min, r_max_val);

        % Color map by cluster ID
        K = max(synID_tmp);
        cmap = hsv(max(K,1));
        cmap = cmap(randi(max(K,1),max(K,1),1),:);
        C = cmap(synID_tmp,:);   % Nevent x 3

        % Draw / update a single scatter (fast)
        if isempty(fig.UserData.scatterHandle) || ~isgraphics(fig.UserData.scatterHandle)
            fig.UserData.scatterHandle = scatter(ax, cent(:,1), cent(:,2), 12, C, 'filled');
        else
            h = fig.UserData.scatterHandle;
            h.XData = cent(:,1);
            h.YData = cent(:,2);
            h.CData = C;
        end

        title(ax, sprintf('Event centroids clustered (r_{max}=%.2f). #clusters=%d', r_max_val, K));
        drawnow;
    end

    function synID = cluster_by_rmin_rmax(cent, r_min, r_max)
        Nevent = size(cent,1);

        % adjacency within r_min
        r2 = r_min^2;
        A = false(Nevent, Nevent);
        for i = 1:Nevent
            dx = cent(i,1) - cent(:,1);
            dy = cent(i,2) - cent(:,2);
            A(i,:) = (dx.^2 + dy.^2) <= r2;
        end

        G = graph(A);
        synID0 = conncomp(G)'; % initial connected components
        clusters = accumarray(synID0, (1:Nevent)', [], @(v){v});

        synID = zeros(Nevent,1);
        synCounter = 0;

        for c = 1:numel(clusters)
            idx = clusters{c};
            pts = cent(idx,:);

            dx = max(pts(:,1)) - min(pts(:,1));
            dy = max(pts(:,2)) - min(pts(:,2));
            bboxDiag = hypot(dx, dy);

            if bboxDiag <= r_max || numel(idx) <= 2
                synCounter = synCounter + 1;
                synID(idx) = synCounter;
                continue
            end

            K = max(2, ceil(bboxDiag / r_max));
            try
                lab = kmeans(pts, K, 'Replicates', 3, 'MaxIter', 200);
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
    end

    function L = fracLimits(v, fLo, fHi, epsPad)
        v = v(isfinite(v));
        v = sort(v(:));
        n = numel(v);
        if n == 0
            L = [0 1];
            return
        end
        iLo = max(1, min(n, floor(fLo*n)));
        iHi = max(1, min(n, ceil(fHi*n)));
        if iHi < iLo, iHi = iLo; end
        L = [v(iLo) v(iHi)];
        if L(1) == L(2)
            L = [L(1) L(1) + epsPad];
        end
    end



end