function trk = flow_track(bp, cfg, ref)
% trk = flow_track(bp, cfg, ref)
%
% Detect red blood cells in the band-passed movie S frame by frame, link them into
% tracklets, and measure how fast and how consistently they travel.
%
% The output that matters downstream is not the trajectories themselves -- it is the
% SPEED per vessel, because the oscillator's frequency band is set by v/d (speed over
% inter-particle spacing). Tracking is used for what it is reliable at: telling you where
% the vessel is and which frequency band to look in.
%
% BRIGHT AND DARK ARE DETECTED SEPARATELY. RBCs are dark against plasma, but the
% band-pass inverts contrast at some depths, so both polarities are run and linked
% independently and each tracklet records which population it came from. Linking them
% together would connect a cell to the plasma gap beside it.
%
% MEMORY IS 0 BY DEFAULT. The Crocker-Grier track() in this repo aborts ("You are not
% going any further") for mem > 0 on dense data. A larger value is attempted and the code
% falls back to 0 automatically, reporting that it did.
%
% THIS STAGE DOES NOT JUDGE. It detects, links, and measures each tracklet's speed,
% length, amplitude and travel consistency -- and keeps all of them. Deciding which are
% trustworthy is the first thing flow_segments does, because that decision defines the
% segments and belongs with them.
%
% INPUT   bp   result of flow_bandpass (uses .S and .varS)
%         cfg  uses cfg.track and cfg.common.fs
%         ref  .pxUm required; .bvBW optional (drawn by the preview only)
% OUTPUT  trk
%   .all    {1 x nTrk} tracklets, columns [x y brightness Rg2 frame]
%   .pol    [nTrk x 1] +1 bright, -1 dark
%   .keep   [nTrk x 1] all true: the quality cuts are applied by flow_segments
%   .cen    [nTrk x 2] centroid, px      .vel [nTrk x 2] velocity, um/s
%   .spd    [nTrk x 1] speed, um/s       .len [nTrk x 1] frames
%   .amp .sig .straight .persist .chanceStr
%   .det    [nDet x 4] every detection before linking: [x y frame polarity]
%   .frames    the frames tracked
%   .speedMap  [nrow x ncol] binned mean speed, for the preview
%   .quiver    [n x 4] x, y, u, v of the binned flow field
%   .nDet      detections per polarity
%
% 2026.08.28, Byung Hun Lee

c  = cfg.track;
fs = cfg.common.fs;
pxUm = ref.pxUm;
[nrow, ncol, nT] = size(bp.S);

frames = c.frames;
if isempty(frames), frames = 1:min(2000, nT); end
frames = frames(frames >= 1 & frames <= nT);
assert(numel(frames) > 10, 'Need more than 10 frames to track; got %d.', numel(frames));

featPx = c.featPx;
if isempty(featPx)                       % derive it from the band the filter achieved
    featPx = 2*round(bp.bandUm(1)/pxUm/2) + 1;
end
featPx = max(3, 2*floor(featPx/2) + 1);  % pkfnd needs an odd size

Strk = bp.S(:,:,frames);
thr  = c.thrSD * std(double(Strk(:)));
switch lower(c.polarity)
    case 'bright', pols = 1;      polNm = {'bright'};
    case 'dark',   pols = -1;     polNm = {'dark'};
    case 'both',   pols = [1 -1]; polNm = {'bright','dark'};
    otherwise, error('cfg.track.polarity must be bright | dark | both');
end
fprintf(['[flow] track: %d frames, threshold %.3g (=%gSD), feature %d px, ' ...
         'max displacement %d px\n'], numel(frames), thr, c.thrSD, featPx, c.maxDisp);

% ---- detect, then link, one polarity at a time ---------------------------------------
trkAll = {};  trkPol = [];  nDet = zeros(1, numel(pols));  nOut = zeros(1, numel(pols));
detAll = zeros(0,4);   % [x y frame polarity] -- every detection, before linking
for pp = 1:numel(pols)
    posT = cell(1, numel(frames));
    for t = 1:numel(frames)
        im = pols(pp) * double(Strk(:,:,t));
        pk = pkfnd(im, thr, featPx);
        if isempty(pk), posT{t} = zeros(0,4); continue, end
        ct = cntrd(im, pk, featPx+2);
        if isempty(ct), posT{t} = zeros(0,4); continue, end
        % cntrd returns an intensity-weighted centroid, sum(x.*I)/sum(I). S is SIGNED
        % (its mean was removed), so a window whose weights nearly cancel sends that
        % ratio anywhere -- measured on real settings, 3.5% of centroids landed outside
        % the frame, some at x = -1355 or y = 4540. They are arithmetic blow-ups, not
        % detections. Linking already ignores them (maxDisp keeps them unlinkable), but
        % they wreck any plot that autoscales to the data, so they are dropped here.
        inFr = ct(:,1) >= 1 & ct(:,1) <= ncol & ct(:,2) >= 1 & ct(:,2) <= nrow;
        nOut(pp) = nOut(pp) + nnz(~inFr);
        ct = ct(inFr, :);
        if isempty(ct), posT{t} = zeros(0,4); continue, end
        posT{t} = ct(:,1:4);                 % [x y brightness Rg^2]
        nDet(pp) = nDet(pp) + size(ct,1);
    end
    if nOut(pp) > 0
        fprintf('       %-6s : %d detections (%.1f per frame), %d centroid blow-ups dropped\n', ...
            polNm{pp}, nDet(pp), nDet(pp)/numel(frames), nOut(pp));
    else
        fprintf('       %-6s : %d detections (%.1f per frame)\n', ...
            polNm{pp}, nDet(pp), nDet(pp)/numel(frames));
    end

    np = cellfun(@(p) size(p,1), posT);
    ppos = zeros(sum(np), 5);
    p = 0;
    for t = 1:numel(posT)
        if np(t) == 0, continue, end
        ppos(p+(1:np(t)), :) = [posT{t}, repmat(t, np(t), 1)];
        p = p + np(t);
    end
    detAll = [detAll; ppos(:,1:2), ppos(:,5), repmat(pols(pp), size(ppos,1), 1)]; %#ok<AGROW>
    if size(ppos,1) < 10
        warning('flow_track:fewDet', 'Too few %s detections to link.', polNm{pp});
        continue
    end

    lk = [];
    for memTry = unique([c.mem 0], 'stable')
        try
            lk = track(ppos, c.maxDisp, struct('mem',memTry,'good',0,'dim',2,'quiet',1));
            if memTry ~= c.mem
                fprintf('       (track() failed at mem=%d, succeeded at mem=%d)\n', c.mem, memTry);
            end
            break
        catch ME
            fprintf('       track() failed at mem=%d: %s\n', memTry, ME.message);
        end
    end
    if isempty(lk)
        warning('flow_track:noLink', ...
            ['track() could not link the %s detections. Reduce maxDisp, raise thrSD so ' ...
             'there are fewer detections, or shorten frames.'], polNm{pp});
        continue
    end
    for idp = unique(lk(:,6))'
        r = lk(lk(:,6)==idp, 1:5);
        if size(r,1) >= c.minLen, trkAll{end+1} = r; trkPol(end+1) = pols(pp); end %#ok<AGROW>
    end
end
nTrk = numel(trkAll);
assert(nTrk > 0, ['No tracklet survived linking. Lower thrSD (more detections), raise ' ...
    'maxDisp, or lower minLen.']);
fprintf('       linked %d tracklets of >= %d frames (%d bright, %d dark)\n', ...
    nTrk, c.minLen, nnz(trkPol > 0), nnz(trkPol < 0));

% ---- per-tracklet geometry ------------------------------------------------------------
vel = nan(nTrk,2);  spd = nan(nTrk,1);  len = nan(nTrk,1);  cen = nan(nTrk,2);
amp = nan(nTrk,1);  sig = nan(nTrk,1);  str = nan(nTrk,1);  per = nan(nTrk,1);
for q = 1:nTrk
    r = trkAll{q};
    amp(q) = median(r(:,3));
    sig(q) = median(sqrt(max(r(:,4),0)/2)) * pxUm;      % Rg^2 = 2 sigma^2
    dxy = (r(end,1:2) - r(1,1:2)) * pxUm;
    dtS = (r(end,5) - r(1,5)) / fs;
    vel(q,:) = dxy / max(dtS, eps);
    spd(q)   = norm(vel(q,:));
    len(q)   = size(r,1);
    cen(q,:) = mean(r(:,1:2), 1);

    st = diff(r(:,1:2), 1, 1);
    sl = vecnorm(st, 2, 2);
    str(q) = norm(r(end,1:2) - r(1,1:2)) / max(sum(sl), eps);
    ok = sl > eps;
    if nnz(ok) >= 2
        u = st(ok,:) ./ sl(ok);
        per(q) = mean(sum(u(1:end-1,:) .* u(2:end,:), 2));
    end
end
chanceStr = sqrt(pi ./ (4*max(len-1,1)));

% ---- how consistently does each tracklet travel? --------------------------------------
% Measured here, JUDGED in flow_segments. Two measures, and the difference matters:
%   straightness = |net displacement| / path length. Its chance level DEPENDS ON LENGTH
%       -- an n-step random walk scores about sqrt(pi/(4n)), so a 6-step tracklet reaches
%       0.36 by luck. Reported next to that chance level, never used as a filter.
%   persistence  = mean cos(angle between consecutive steps). 0 for a random walk
%       REGARDLESS of length, 1 for straight travel. Length-fair, so this is what
%       flow_segments filters on.
fprintf(['       speed median %.0f um/s (IQR %.0f-%.0f) | length median %d frames\n' ...
         '       persistence %+.3f (0 = random walk) | straightness %.3f vs chance %.3f\n'], ...
    median(spd), prctile(spd,25), prctile(spd,75), round(median(len)), ...
    median(per,'omitnan'), median(str,'omitnan'), median(chanceStr));
fprintf('       (the quality cuts are in step 3, flow_segments)\n');
keep = true(nTrk,1);        % nothing is rejected here

% ---- binned flow field, from the kept tracklets only -----------------------------------
xb = 1:c.binPx:ncol;  yb = 1:c.binPx:nrow;
Ub = zeros(numel(yb), numel(xb));  Vb = Ub;  Nb = Ub;  Sb = Ub;
for q = find(keep)'
    bi = max(1, min(numel(yb), ceil(cen(q,2)/c.binPx)));
    bj = max(1, min(numel(xb), ceil(cen(q,1)/c.binPx)));
    Ub(bi,bj) = Ub(bi,bj) + vel(q,1);
    Vb(bi,bj) = Vb(bi,bj) + vel(q,2);
    Sb(bi,bj) = Sb(bi,bj) + spd(q);
    Nb(bi,bj) = Nb(bi,bj) + 1;
end
Ub = Ub ./ max(Nb,1);  Vb = Vb ./ max(Nb,1);  Sb = Sb ./ max(Nb,1);
okB = Nb >= 2;                      % a bin needs two tracklets before it means anything

% paint the binned speed back onto the full frame, for the preview
speedMap = zeros(nrow, ncol);
for bi = 1:numel(yb)
    for bj = 1:numel(xb)
        if ~okB(bi,bj), continue, end
        ri = yb(bi):min(yb(bi)+c.binPx-1, nrow);
        ci = xb(bj):min(xb(bj)+c.binPx-1, ncol);
        speedMap(ri,ci) = Sb(bi,bj);
    end
end
[Xb, Yb] = meshgrid(xb + c.binPx/2, yb + c.binPx/2);
qs = c.binPx / max(prctile(spd(keep), 90), eps);       % arrows scaled to fit a bin
quiv = [Xb(okB), Yb(okB), Ub(okB)*qs, Vb(okB)*qs];

trk = struct('all', {trkAll}, 'pol', trkPol(:), 'keep', keep, 'cen', cen, ...
    'vel', vel, 'spd', spd, 'len', len, 'amp', amp, 'sig', sig, ...
    'straight', str, 'persist', per, 'chanceStr', chanceStr, ...
    'frames', frames, 'featPx', featPx, 'thr', thr, 'nDet', nDet, 'nOutOfFrame', nOut, 'det', detAll, ...
    'speedMap', speedMap, 'quiver', quiv, 'quiverScale', qs);
end
