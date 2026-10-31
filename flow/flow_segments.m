function seg = flow_segments(trk, sz, cfg, ref)
% seg = flow_segments(trk, sz, cfg, ref)
%
% Turn the tracklets into VESSEL SEGMENTS: where were RBCs actually found, which
% connected patches of that are big enough to be a vessel, and which tracklet belongs to
% which patch.
%
% WHAT THIS REPLACES. Cell (7e) of the script ran 253 lines, but only about fifteen of
% them fed anything downstream. Its headline product -- P, a movie rendered from the
% detections, and the r2map asking how much of S that rendering explains -- was a
% NEGATIVE result: the rendering is zero wherever detection failed, and as a regressor it
% reached R^2 = 0.30 against 0.82 for the phase oscillator that replaced it. Cell (7f)
% consumed only densBW, detXY, detT and trkUse; detA, detS, P and r2map appeared
% downstream solely inside a comment and an assert message.
%
% That comparison lives in cell (7e) of the original script, which is kept for the record.
% It is NOT reimplemented here: this stage produces the segment mask and nothing else.
%
% NO NEW DETECTION IS RUN. detXY/detT are a flattening of the tracklets flow_track
% already built -- the script's apparent second detection pass was only that.
%
% INPUT
%   trk   result of flow_track: .all {1 x nTrk} tracklets [x y bright Rg2 frame]
%                               .pol .cen .spd, and .len/.amp/.persist if available
%   sz    [nrow ncol]
%   cfg   full config; uses cfg.segments
%   ref   optional struct of reference masks for the printout only:
%           .bvBW  logical vessel mask   .pxUm  um per pixel
%
% OUTPUT  seg
%   .px        {1 x nSeg}  pixel indices of each USED segment
%   .useIdx    [1 x nSeg]  index of each used segment in the full connected-component list
%   .of        [nrow x ncol] segment label per pixel (0 = none)
%   .densBW    logical mask the segments came from
%   .densMap   raw detection-count image
%   .detXY     [nDet x 2] detection positions   .detT [nDet x 1] frame
%   .detSeg    [nDet x 1] which segment each detection fell in (0 = none)
%   .trkSeg    [nTrk x 1] which segment each tracklet's centroid fell in
%   .trkUse    indices of the tracklets that passed the sort
%   .keep      [nTrk x 1] logical, the sort result
%   .cut       struct of the three individual cuts, for the preview
%   .stat      struct with len/amp/persist per tracklet
%   .nTrkSeg   [1 x nSeg] tracklets per used segment
%   .segSpeed  [1 x nSeg] median tracklet speed in each used segment, um/s. This is the
%              v that flow_regressors turns into a frequency band via v/d.
%   .segNet    [1 x nSeg] median net travel per segment, um
%   .trkSpd    [nTrk x 1] speed of every tracklet, um/s
%
% 2026.08.28, Byung Hun Lee

if nargin < 4, ref = struct; end
c = cfg.segments;
nrow = sz(1);  ncol = sz(2);

% ---- STEP 1: SORT THE TRACKLETS -----------------------------------------------------
% Nothing downstream can recover from a bad choice here: only the surviving tracklets
% contribute detections, so these three cuts decide where the density map -- and
% therefore every segment -- ends up.
%
% NET TRAVEL FIRST -- how far the cell actually got, first point to last. It is the most
% direct statement of "this thing moved", it is in physical units so it can be set from
% what a red cell is expected to do, and it fails a tracklet that wandered in place
% however long that tracklet is. Note it is the net DISPLACEMENT, not the path length: a
% cell jittering back and forth covers a long path while going nowhere.
%
% Duration, amplitude and persistence are secondary. Persistence rather than
% straightness, because straightness scores 0.36 by luck on a 6-step random walk while
% persistence is 0 for a random walk of ANY length.
%
% flow_track measures len/amp/persist and keeps everything; the numbers are recomputed
% here from trk.all when a caller supplies tracklets without them.
nTrkAll = numel(trk.all);
assert(nTrkAll > 0, 'flow_track returned no tracklets.');
pxUm = 1;  if isfield(ref,'pxUm') && ~isempty(ref.pxUm), pxUm = ref.pxUm; end
fsHz = 1000; if isfield(ref,'fs') && ~isempty(ref.fs), fsHz = ref.fs; end
[tLen, tAmp, tPer, tNet] = trackStats(trk, nTrkAll, pxUm);

cutNet = tNet >= c.minNetUm;
cutLen = tLen >= c.minLenKeep;
cutAmp = tAmp >= prctile(tAmp, c.minAmpPrc);
cutPer = tPer >= c.minPersist | isnan(tPer);
keep   = cutNet & cutLen & cutAmp & cutPer;
fprintf('[flow] segments, tracklet sort (of %d):\n', nTrkAll);
fprintf('       net travel  >= %4.1f um   : %4d pass  <- primary\n', c.minNetUm, nnz(cutNet));
fprintf('       duration    >= %2d frames : %4d pass\n', c.minLenKeep, nnz(cutLen));
fprintf('       |amp|       >= p%-4g      : %4d pass\n', c.minAmpPrc, nnz(cutAmp));
fprintf('       persistence >= %+.2f      : %4d pass\n', c.minPersist, nnz(cutPer));
fprintf('       ALL FOUR                 : %4d kept  (net travel median %.1f um)\n', ...
    nnz(keep), median(tNet(keep)));
assert(any(keep), ['Every tracklet was rejected. The cut that usually does this is ' ...
    'cfg.segments.minNetUm (currently %.1f um) or minLenKeep (%d frames).'], ...
    c.minNetUm, c.minLenKeep);
if nnz(keep) < 10
    warning('flow_segments:fewKept', ...
        ['Only %d tracklets survived the sort; the density map will be sparse and the ' ...
         'segments unreliable. Lower minLenKeep or cfg.track.thrSD.'], nnz(keep));
end
trkUse = find(keep)';

% ---- flatten them into detections ---------------------------------------------------
nDet = sum(cellfun(@(r) size(r,1), trk.all(trkUse)));
detXY = zeros(nDet,2);  detT = zeros(nDet,1);  detA = zeros(nDet,1);  detS = zeros(nDet,1);
p = 0;
for q = trkUse
    r  = trk.all{q};  n = size(r,1);
    ii = p + (1:n);
    detXY(ii,:) = r(:,1:2);
    detA(ii)    = trk.pol(q) * r(:,3);          % signed integrated brightness
    detS(ii)    = sqrt(max(r(:,4),0)/2);        % sigma in pixels (Rg^2 = 2 sigma^2)
    detT(ii)    = r(:,5);
    p = p + n;
end

% ---- detection-density mask ---------------------------------------------------------
lin = sub2ind([nrow ncol], min(max(round(detXY(:,2)),1),nrow), ...
                           min(max(round(detXY(:,1)),1),ncol));
densMap = accumarray(lin, 1, [nrow*ncol 1]);
densMap = reshape(densMap, nrow, ncol);
densSm  = imgaussfilt(densMap, c.densSm);
pos     = densSm(densSm > 0);
assert(~isempty(pos), 'The density map is empty: no detection landed anywhere.');
densBW  = densSm > prctile(pos, c.densPrc);

% ---- connected components -> candidate segments -------------------------------------
CC    = bwconncomp(densBW);
segPx = CC.PixelIdxList;
segPx = segPx(cellfun(@numel, segPx) >= c.minPx);
nSeg  = numel(segPx);
assert(nSeg > 0, ['No connected patch reached minPx = %d. Lower cfg.segments.densPrc ' ...
    'so the mask is bigger, or minPx so smaller patches qualify.'], c.minPx);

segOf = zeros(nrow, ncol);
for s = 1:nSeg, segOf(segPx{s}) = s; end

detSeg = segOf(lin);
trkSeg = zeros(numel(trk.all), 1);
for q = trkUse
    cl = sub2ind([nrow ncol], min(max(round(trk.cen(q,2)),1),nrow), ...
                              min(max(round(trk.cen(q,1)),1),ncol));
    trkSeg(q) = segOf(cl);
end

% ---- a segment needs enough tracklets to give a trustworthy speed --------------------
nTrkSeg = arrayfun(@(s) nnz(trkSeg == s), 1:nSeg);
useIdx  = find(nTrkSeg >= c.minTrk);
assert(~isempty(useIdx), ['No segment reached minTrk = %d (best %d). Either the ' ...
    'tracking found too little, or the density mask is cutting vessels into pieces.'], ...
    c.minTrk, max(nTrkSeg));

% ---- per-segment speed ---------------------------------------------------------------
% Stored here, not recomputed later: flow_regressors sets each segment's frequency band
% from v/d, so this median speed IS the v in that expression. Taking the median rather
% than the mean matters -- one mislinked tracklet can carry an absurd velocity, and the
% band it would produce is nonsense.
if isfield(trk,'spd') && numel(trk.spd) == nTrkAll
    tSpd = trk.spd(:);
else
    tSpd = nan(nTrkAll,1);
    for q = 1:nTrkAll
        r = trk.all{q};
        dt = (r(end,5) - r(1,5)) / max(fsHz, eps);
        tSpd(q) = norm(r(end,1:2) - r(1,1:2)) * pxUm / max(dt, eps);
    end
end
segSpeed = arrayfun(@(s) median(tSpd(trkSeg == s), 'omitnan'), useIdx);
segNet   = arrayfun(@(s) median(tNet(trkSeg == s), 'omitnan'), useIdx);

seg = struct('keep', keep, 'cut', struct('net',cutNet,'len',cutLen,'amp',cutAmp,'per',cutPer), ...
    'segSpeed', segSpeed, 'segNet', segNet, 'trkSpd', tSpd, ...
    'stat', struct('net',tNet,'len',tLen,'amp',tAmp,'persist',tPer), ...
    'px', {segPx(useIdx)}, 'useIdx', useIdx, 'of', segOf, ...
    'densBW', densBW, 'densMap', densMap, 'detXY', detXY, 'detT', detT, ...
    'detA', detA, 'detS', detS, 'detSeg', detSeg, 'trkSeg', trkSeg, ...
    'trkUse', trkUse, 'nTrkSeg', nTrkSeg(useIdx), 'nSegAll', nSeg);

fprintf('[flow] segments: %d detections from %d tracklets -> %d px mask (p%g)\n', ...
    nDet, numel(trkUse), nnz(densBW), c.densPrc);
fprintf('       %d connected patches >= %d px, %d with >= %d tracklets\n', ...
    nSeg, c.minPx, numel(useIdx), c.minTrk);
if isfield(ref,'bvBW') && ~isempty(ref.bvBW)
    fprintf('       %.0f%% of the mask falls inside bvMask\n', 100*mean(ref.bvBW(densBW)));
end

end

function [tLen, tAmp, tPer, tNet] = trackStats(trk, n, pxUm)
% length, amplitude and persistence per tracklet -- from flow_track when it supplied
% them, recomputed from the raw tracklets otherwise, so a caller can pass in tracklets
% from anywhere (the script's trkAll, for instance) and still be sorted the same way.
if isfield(trk,'len') && numel(trk.len) == n, tLen = trk.len(:);
else, tLen = cellfun(@(r) size(r,1), trk.all(:));
end
if isfield(trk,'amp') && numel(trk.amp) == n, tAmp = trk.amp(:);
else, tAmp = cellfun(@(r) median(r(:,3)), trk.all(:));
end
if isfield(trk,'persist') && numel(trk.persist) == n
    tPer = trk.persist(:);
else
    tPer = nan(n,1);
    for q = 1:n
        st = diff(trk.all{q}(:,1:2), 1, 1);
        sl = vecnorm(st, 2, 2);
        ok = sl > eps;
        if nnz(ok) >= 2
            u = st(ok,:) ./ sl(ok);
            tPer(q) = mean(sum(u(1:end-1,:) .* u(2:end,:), 2));
        end
    end
end
% net displacement, first point to last, in um. Deliberately NOT the path length: a
% cell that jitters back and forth covers a long path while going nowhere.
tNet = cellfun(@(r) norm(r(end,1:2) - r(1,1:2)), trk.all(:)) * pxUm;
end
