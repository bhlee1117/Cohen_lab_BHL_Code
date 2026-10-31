function reg = flow_regressors(bp, trk, seg, data, cfg)
% reg = flow_regressors(bp, trk, seg, data, cfg)
%
% Build the vascular design matrix: one phase-continuous oscillator per vessel segment.
%
% =====================================================================================
% THE IDEA
%
% Red cells pass a point in a capillary quasi-periodically, so the local signal is
%
%       s(t) = A(t) * cos(phi(t)),      dphi/dt = 2*pi*f(t)
%
% and a regressor written that way keeps running through the gaps where detection
% failed. That is the whole reason this replaced the obvious alternative: rendering the
% tracked particles back into a movie gives a regressor that is ZERO wherever a cell was
% missed, punching holes exactly where the vascular signal is still present. On a
% synthetic segment with 34% dropout the rendering reached R^2 = 0.30; the phase
% oscillator reached 0.82.
%
% HOW *NOT* TO GET f(t). The obvious estimate is the particle passage rate, detections
% per second. It fails: every missed detection lowers the count, and phi is the running
% INTEGRAL of f, so the error accumulates until the oscillator is out of phase with the
% signal. On the same synthetic test a true 40 Hz came out as 20 Hz and explained
% R^2 = 0.17 -- worse than the rendering.
%
% WHAT WORKS. Take the phase and amplitude from the MEASURED segment signal (band-pass +
% Hilbert), and use the tracking only for what it is reliable at: WHERE the vessel is,
% and WHICH BAND to look in, via v/d -- speed over inter-cell spacing, neither of which
% depends on how many other cells were detected. The band matters: centred right with a
% +/-60% width it works, while a narrow band centred on the biased rate scores 0.06.
%
% =====================================================================================
% THE STEPS, per segment
%
%  1. BAND FROM THE FLOW.  f = v/d, with v the median tracklet speed in that segment
%     (flow_segments measured it) and d the median nearest-neighbour distance between
%     SIMULTANEOUS detections. Search band = f * [1/3, 3], then clamped to the floor.
%
%  2. DECOMPOSE THE SEGMENT.  A segment longer than the inter-cell spacing holds an
%     alternating dark/bright train, so its spatial MEAN cancels -- on a synthetic
%     vessel 6x longer than the spacing the mean suppressed the envelope 68-fold. So the
%     segment is decomposed instead, by COMPLEX PCA: band-limit each pixel, take its
%     analytic signal, then SVD. A travelling wave is rank 2 in the real domain (two
%     modes 90 deg apart, its variance split between them) and rank 1 in the complex
%     one -- measured on a synthetic 60 Hz wave, real PCA gave 49%/49% across two modes
%     while complex PCA put 96% in one and recovered the phase gradient at corr 1.000.
%     |U| is then the envelope along the vessel and angle(U) is the propagation.
%
%  3. KEEP THE OSCILLATORY MODES.  get_oscScore = spectral concentration x
%     autocorrelation side lobe. Both are needed: concentration alone scores white noise
%     0.29, the side lobe alone scores pure drift 0.999, the product separates them
%     (sine 0.97, drift 0.04, white 0.02). Modes are ranked by score x variance, so an
%     oscillation nobody can see does not outrank a strong one.
%
%  4. PHASE.  PSD peak of the leading mode inside the search band -> fPeak; passband
%     fPeak*(1 +/- bandFrac); band-pass, Hilbert, and emit BOTH cos and sin. Being in
%     quadrature, a per-pixel fit on the pair matches any phase offset between the
%     segment and that pixel, so no alignment step is needed.
%
%  5. THE FLOOR.  Every frequency here is clamped below by cfg.regressors.fMin.
%     Dendritic Vm lives under it; RBC passage (v ~ 0.5-2 mm/s, d ~ 5-15 um) is
%     40-400 Hz. The floor clamps the search, the passband AFTER the +/-bandFrac
%     widening, and the fallback -- all three, because a legal fPeak of 50 Hz would
%     otherwise open the band down to 20 Hz. A segment whose peak lands ON the floor had
%     its band imposed rather than measured and is dropped.
%
%  6. TWO TESTS ON EVERY ROW, AND THEY ARE NOT THE SAME KIND OF TEST.  oscScore grades
%     how much the row looks like an oscillation IN TIME, and scoreMin is a threshold on
%     that grade. The voltage map asks WHERE it lives, and that one is a VETO: each row's
%     spatial map is correlated against a voltage reference, and |r| >= voltCut drops the
%     row outright, whatever it scored. A periodic thing sitting exactly on the dendrite
%     is far more likely to be voltage the design has latched onto than blood, and there
%     is no |r| at which eating some voltage becomes an acceptable trade -- so there is
%     nothing for a graded penalty to express. Rows killed this way are listed in
%     .vetoed, since they are not in the output to be looked at afterwards.
% =====================================================================================
%
% INPUT   bp   flow_bandpass result (S, mov_nv)
%         trk  flow_track result
%         seg  flow_segments result (px, segSpeed, detXY/detT/detSeg)
%         data .mov_res .badfr .spike .neuronMap .ftprnt .pxUm
%         cfg  uses cfg.regressors and cfg.common
% OUTPUT  reg
%   .design   [nReg x nT] the design matrix, z-scored rows
%   .names    {1 x nReg}  what each row is
%   .segOf    [nReg x 1]  which entry of .useSeg each row came from
%   .maps     [nPix x nReg] covariance of each row with mov_res: its spatial footprint
%   .score .oscScore .rVolt .rVoltAlt .f0   per row (.score IS .oscScore: the voltage
%             test no longer scales it, it vetoes instead)
%   .vetoed   what the voltage veto removed: .name .rVolt .rVoltAlt .oscScore .cut
%   .varExpl  [nReg x 1] % of the movie variance this row explains ON ITS OWN
%   .varFit   [nrow x ncol] variance the WHOLE design removes at each pixel (joint fit)
%   .dropFrac [nrow x ncol] that, as a fraction of the movie's own variance
%   .varMov   [nrow x ncol] the movie's variance, for reference
%   .eigImg   [nPix x nReg] complex spatial mode behind each row (NaN off-segment)
%   .modeOf   [nReg x 1] which mode of its segment each row came from
%   .useSeg   segments that survived            .frames  frames the design spans
%   .seg      per-segment struct array: band, fPeak, fTrk, spacing, mode, U, srcTrace
%   .voltRef  which reference the veto used, and the map itself
%
% 2026.08.28, Byung Hun Lee

c   = cfg.regressors;
fs  = cfg.common.fs;
frames = trk.frames;
nF  = numel(frames);
[nrow, ncol, ~] = size(bp.S);
npx = nrow*ncol;
pxUm = data.pxUm;
fFloor = max(c.fMin, 1);

% nanPSD parallelises over its rows; it is called here with ONE row per segment, so the
% pool is pure overhead (46 s vs 2.6 s measured). Restored when this function returns.
poolGuard = flow_noPool();  %#ok<NASGU>

gF = ~data.badfr(frames);                 % no spike, no blue: the frames every fit uses
assert(nnz(gF) > 50, 'Only %d good frames in the tracked window.', nnz(gF));

nSeg = numel(seg.px);
Sv   = reshape(bp.S, [], size(bp.S,3));   % a view; only segment pixels are ever copied
Mnv  = reshape(bp.mov_nv, [], size(bp.mov_nv,3));

% ---- per segment: band, decomposition, phase ----------------------------------------
oscCos = zeros(nSeg, nF);  oscSin = oscCos;  oscY = oscCos;  srcNv = oscCos;
sInfo(nSeg) = struct('band',[], 'fPeak',[], 'fTrk',[], 'spacing',[], 'mode',[], ...
                     'U',[], 'pcIdx',[], 'nPC',[], 'pcScore',[], 'pcAll',[], 'srch',[], ...
                     'srcFull',[]);   % the S-derived trace the PSD and the phase came from
pinned = false(1,nSeg);
fprintf('[flow] regressors: %d segments, %s, floor %g Hz\n', nSeg, ...
    ternary(c.cpca,'complex PCA','PCA'), c.fMin);

for k = 1:nSeg
    pxIdx = seg.px{k};
    sgi   = seg.useIdx(k);

    % --- 1. the band the flow implies -------------------------------------------------
    vMed = seg.segSpeed(k);
    dsel = seg.detSeg == sgi;
    dd = [];
    tq = seg.detT(dsel);  xy = seg.detXY(dsel,:);
    for t = unique(round(tq))'
        p = xy(round(tq) == t, :);
        if size(p,1) < 2, continue, end
        D = squareform(pdist(p)) + diag(inf(size(p,1),1));
        dd = [dd; min(D,[],2)]; %#ok<AGROW>
    end
    spacing = median(dd) * pxUm;
    fTrk = vMed / max(spacing, eps);

    fLim = [fFloor, 0.45*fs];
    if ~isfinite(fTrk) || fTrk < fLim(1) || fTrk > fLim(2)
        fprintf(['   seg %3d: v/d = %.3g Hz is outside %g-%g Hz (spacing %.2f um, ' ...
                 'speed %.0f um/s) -> searching the full band\n'], ...
            sgi, fTrk, fLim(1), fLim(2), spacing, vMed);
        srch = fLim;
    else
        srch = [fTrk/c.search(2), fTrk*c.search(2)];
    end
    srch = [max(srch(1), fFloor), min(srch(2), fs/2 - 1)];
    if srch(2) <= srch(1)*1.05, srch = [fFloor, fs/2 - 1]; end

    % --- 2/3. decompose and keep the oscillatory modes ---------------------------------
    Y = -double(Sv(pxIdx, frames));
    Y = Y - mean(Y, 2);
    if numel(pxIdx) >= c.minPxPCA
        nScan = max(2, ceil(numel(pxIdx)/c.binPix));
        nScan = min(nScan, min(size(Y,1), nnz(gF)) - 1);
        if c.cpca
            [U, pcTr, varPC] = get_cpca(Y, fs, srch, nScan, gF);
            scPC = get_oscScore(real(pcTr), fs, srch);
        else
            [U, Ss, ~] = svds(Y(:, gF), nScan);
            pcTr  = U' * Y;
            varPC = diag(Ss).^2;  varPC = varPC / max(sum(varPC), eps);
            scPC  = get_oscScore(pcTr, fs, srch);
        end
        keepPC = find(scPC >= c.scoreMin);
        [~, ord] = sort(scPC(keepPC) .* varPC(keepPC), 'descend');   % oscillatory AND strong
        keepPC = keepPC(ord);
        if numel(keepPC) > c.maxPC, keepPC = keepPC(1:c.maxPC); end
        if isempty(keepPC)
            warning('flow_regressors:noPC', ...
                ['Segment %d: no mode reached scoreMin = %.2f (best %.2f); falling back ' ...
                 'to the segment mean, which may be phase-cancelled.'], ...
                sgi, c.scoreMin, max(scPC));
            y = mean(Y,1);  md = 2;  U = [];  keepPC = [];  pcKeep = [];
        else
            pcKeep = pcTr(keepPC,:);
            y = real(pcKeep(1,:));  md = 1;  U = U(:,keepPC);
        end
        fprintf('   seg %3d: %4d px  %d modes -> %d kept (score %s | var %s%%)\n', ...
            sgi, numel(pxIdx), numel(scPC), numel(keepPC), ...
            mat2str(round(scPC(keepPC)',2)), mat2str(round(100*varPC(keepPC)')));
    else
        y = mean(Y,1);  md = 3;  U = [];  keepPC = [];  pcKeep = [];  scPC = [];
        fprintf('   seg %3d: %4d px  own trace (< minPxPCA = %d)\n', ...
            sgi, numel(pxIdx), c.minPxPCA);
    end

    % --- 4/5. the peak, the passband, the phase ---------------------------------------
    yg = y;  yg(~gF) = NaN;
    [Py, Fy] = nanPSD(yg, fs, min(1024, 2^floor(log2(nF/2))));
    Py = Py(:);  Fy = Fy(:);
    in = Fy >= srch(1) & Fy <= srch(2);
    if ~any(in & isfinite(Py)), in = Fy > fFloor & Fy < fs/2-1; end
    Pm = Py;  Pm(~in | ~isfinite(Pm)) = -inf;
    [~, ip] = max(Pm);  fPeak = Fy(ip);
    pinned(k) = c.fMin > 0 && fPeak <= fFloor*1.05;
    if pinned(k)
        fprintf(['   seg %3d: PSD peak pinned at the %g Hz floor -- its real peak is ' ...
                 'below fMin%s\n'], sgi, fFloor, ternary(c.dropPinned,' -> DROPPED',''));
    end

    bnd = fPeak * [1-c.bandFrac, 1+c.bandFrac];
    bnd = [max(bnd(1), fFloor), min(bnd(2), fs/2-1)];      % floor AFTER the widening
    if ~all(isfinite(bnd)) || bnd(2) <= bnd(1)*1.05, bnd = srch; end

    yb = get_bandpass(y, fs, bnd);  yb(~isfinite(yb)) = 0;
    z  = hilbert(yb);
    oscCos(k,:) = abs(z).*cos(angle(z));
    oscSin(k,:) = abs(z).*sin(angle(z));
    oscY(k,:)   = y;

    % the browser's "before any filtering" trace: the SAME spatial weights, applied to
    % mov_nv. y above is already through the spatial band-pass, so showing it as the
    % source would hide what that filter did.
    Ynv = -double(Mnv(pxIdx, frames));
    Ynv = Ynv - mean(Ynv, 2);
    if md == 1 && ~isempty(U), srcNv(k,:) = real(U(:,1)' * Ynv);
    else,                      srcNv(k,:) = mean(Ynv, 1);
    end

    sInfo(k) = struct('band',bnd, 'fPeak',fPeak, 'fTrk',fTrk, 'spacing',spacing, ...
        'mode',md, 'U',U, 'pcIdx',keepPC, 'nPC',numel(scPC), ...
        'pcScore',scPC(keepPC), 'pcAll',pcKeep, 'srch',srch, 'srcFull',y);
end

% ---- drop the segments pinned on the floor ------------------------------------------
% Their true spectral maximum is below fMin, so the band was imposed rather than found,
% and on a dendrite what the lowest legal band contains is voltage.
if c.dropPinned && any(pinned)
    kp = ~pinned;
    fprintf('[flow] dropping %d of %d segments pinned at the %g Hz floor: %s\n', ...
        nnz(pinned), nSeg, fFloor, mat2str(seg.useIdx(pinned)));
    oscCos = oscCos(kp,:);  oscSin = oscSin(kp,:);  oscY = oscY(kp,:);
    srcNv  = srcNv(kp,:);   sInfo  = sInfo(kp);
    useSeg = seg.useIdx(kp);  segPxU = seg.px(kp);
    assert(~isempty(useSeg), ['Every segment pinned at fMin = %g Hz: no vascular ' ...
        'oscillation above the floor was found anywhere.'], c.fMin);
elseif any(pinned)
    warning('flow_regressors:pinned', ...
        '%d segment(s) pin at the %g Hz floor and dropPinned is false.', nnz(pinned), fFloor);
    useSeg = seg.useIdx;  segPxU = seg.px;
else
    useSeg = seg.useIdx;  segPxU = seg.px;
end
nS = numel(useSeg);

% ---- how much of the segment signal does its own oscillator explain? -----------------
r2 = nan(1,nS);
for k = 1:nS
    yv = oscY(k,gF)';  yv = yv - mean(yv);
    X1 = [oscCos(k,gF); oscSin(k,gF); ones(1,nnz(gF))]';
    r2(k) = 1 - var(yv - X1*(X1\yv)) / max(var(yv), eps);
end
fprintf('[flow] oscillator explains its own segment: median R^2 = %.3f (%s)\n', ...
    median(r2), mat2str(round(r2,2)));

% ---- assemble the design -------------------------------------------------------------
zrow = @(V) (V - mean(V,2)) ./ max(std(V,0,2), eps);
srcNm = cell(1,nS);
for k = 1:nS
    if sInfo(k).mode == 1
        tag = ternary(c.cpca,'cpcaPC','svdPC');
        srcNm{k} = sprintf('seg%d %s%d (%d of %d)', useSeg(k), tag, ...
            sInfo(k).pcIdx(1), numel(sInfo(k).pcIdx), sInfo(k).nPC);
    elseif sInfo(k).mode == 2
        srcNm{k} = sprintf('seg%d mean (no mode of %d)', useSeg(k), sInfo(k).nPC);
    else
        srcNm{k} = sprintf('seg%d own trace (small)', useSeg(k));
    end
end
design = [zrow(oscCos); zrow(oscSin)];
names  = [cellfun(@(v)[v ' | cos'], srcNm, 'uni',0), ...
          cellfun(@(v)[v ' | sin'], srcNm, 'uni',0)];
segOf  = [1:nS, 1:nS];
modeOf = [ones(1,nS), ones(1,nS)];    % which spatial mode each row's image comes from
hasQuad = true(1, 2*nS);              % cos and sin ARE each other's 90-degree partner

% Second and later modes are genuinely different spatial patterns, so they go in too.
% Under complex PCA each is its own quadrature pair, hence real AND imag.
for k = 1:nS
    P = sInfo(k).pcAll;
    if isempty(P) || size(P,1) < 2, continue, end
    ex = P(2:end,:);  q = 2:size(P,1);
    if c.cpca
        design = [design; zrow(real(ex)); zrow(imag(ex))];                      %#ok<AGROW>
        names  = [names, ...
            arrayfun(@(j) sprintf('%s | extra%d re', srcNm{k}, sInfo(k).pcIdx(j)), q, 'uni',0), ...
            arrayfun(@(j) sprintf('%s | extra%d im', srcNm{k}, sInfo(k).pcIdx(j)), q, 'uni',0)]; %#ok<AGROW>
        segOf  = [segOf, repmat(k,1,2*size(ex,1))];                             %#ok<AGROW>
        modeOf = [modeOf, q, q];                                                %#ok<AGROW>
        hasQuad = [hasQuad, true(1,2*size(ex,1))];   % re and im are a quadrature pair
    else
        design = [design; zrow(ex)];                                            %#ok<AGROW>
        names  = [names, arrayfun(@(j) sprintf('%s | extra%d', srcNm{k}, ...
            sInfo(k).pcIdx(j)), q, 'uni',0)];                                   %#ok<AGROW>
        segOf  = [segOf, repmat(k,1,size(ex,1))];                               %#ok<AGROW>
        modeOf = [modeOf, q];                                                   %#ok<AGROW>
        hasQuad = [hasQuad, false(1,size(ex,1))];    % a real mode has no partner
    end
end
% derivatives: mostly redundant, since d/dt(A cos) is the sin term weighted by the
% instantaneous frequency (corr -1.000 at constant frequency, -0.964 at 40% FM)
dC = oscCos; dS = oscSin;
for dv = 1:c.useDeriv
    dC = gradient(dC')';  dS = gradient(dS')';
    design = [design; zrow(dC); zrow(dS)];                                      %#ok<AGROW>
    names  = [names, cellfun(@(v)sprintf('%s | d%d cos',v,dv), srcNm,'uni',0), ...
                     cellfun(@(v)sprintf('%s | d%d sin',v,dv), srcNm,'uni',0)]; %#ok<AGROW>
    segOf  = [segOf, 1:nS, 1:nS];                                               %#ok<AGROW>
    modeOf = [modeOf, ones(1,nS), ones(1,nS)];                                  %#ok<AGROW>
    hasQuad = [hasQuad, true(1,2*nS)];      % d/dt cos and d/dt sin pair up too
end

% ---- optional: a 90-degree twin for any row that lacks one --------------------------
if c.addQuad
    need = find(~hasQuad);
    if isempty(need)
        fprintf(['[flow] addQuad: every row already has a 90-degree partner (cos/sin, and ' ...
                 'real/imag\n       under complex PCA), so nothing was added. A quadrature ' ...
                 'pair spans every\n       phase, so if vessel signal survives, phase ' ...
                 'coverage is not the cause.\n']);
    else
        Q = imag(hilbert(design(need,:).').');      % hilbert works down columns
        design = [design; zrow(Q)];                                             %#ok<AGROW>
        names  = [names, cellfun(@(v)[v ' | quad'], names(need), 'uni', 0)];    %#ok<AGROW>
        segOf  = [segOf, segOf(need)];                                          %#ok<AGROW>
        modeOf = [modeOf, modeOf(need)];                                        %#ok<AGROW>
        fprintf('[flow] addQuad: added %d 90-degree twins for rows with no partner\n', numel(need));
    end
end

% ---- each row's spatial footprint, from mov_res --------------------------------------
% Covariance rather than a fitted beta: the rows are correlated with each other, and a
% beta map splits shared structure between them in a way that is not interpretable per
% row. Covariance is marginal -- "where in the movie does THIS waveform appear".
% It is taken from mov_res, which still HAS the voltage in it: a map made from mov_nv
% could not show a regressor tracking voltage, because that voltage is already gone.
Mv = single(reshape(data.mov_res(:,:,frames), npx, nF));
Mv = Mv(:,gF);  Mv = Mv - mean(Mv,2);
Rz = design(:,gF) - mean(design(:,gF),2);
maps = double(Mv * single(Rz)') / max(nnz(gF)-1,1);
clear Mv Rz

% ---- how much of the movie does each row explain, ON ITS OWN? -------------------------
% Marginal, not a joint fit: the rows are correlated with each other, so a joint beta
% splits shared structure between them and cannot be read per row. Each design row is
% z-scored, so its univariate coefficient at a pixel IS the covariance already computed
% above, and the variance it explains there is that covariance squared. Summing over
% pixels and dividing by the movie's own variance gives a fraction that is comparable
% between rows and is what the regressors are ordered by.
% The denominator is the movie's variance INSIDE THE FOOTPRINT MASK. Over the whole
% frame it would be dominated by background pixels that carry no signal, and every
% regressor would read as explaining a fraction of a percent.
aM   = data.anaMask(:);
vPx  = double(var(single(reshape(data.mov_res(:,:,frames(gF)), npx, nnz(gF))), 0, 2));
vTot = sum(vPx(aM));
varExpl = 100 * sum(maps(aM,:).^2, 1)' / max(vTot, eps);

% ---- the spatial mode behind each row, painted on the frame --------------------------
% Complex under CPCA: |U| is the envelope along the vessel, angle(U) the propagation.
% Kept complex here so the browser can show phase-by-amplitude rather than throwing the
% propagation away.
eigImg = nan(npx, size(design,1));
for r = 1:size(design,1)
    k = segOf(r);  j = modeOf(r);
    U = sInfo(k).U;
    if isempty(U) || size(U,2) < j, continue, end
    u = U(:,j);  u = u / max(abs(u));
    eigImg(segPxU{k}, r) = u;
end

% ---- score in time, then penalise in space -------------------------------------------
[oscScore, f0] = get_oscScore(design, fs, [1 fs/2-1]);
[vRef, vName, vAlt] = voltageRef(data, cfg, fs);
% WHERE the correlation is taken: data.anaMask, the same mask every variance figure in
% the pipeline uses. The question is "does this regressor's map look like the neuron",
% so it is asked over the pixels the analysis is about. Widening it to background dilutes
% the comparison with pixels where the voltage reference is ~0 by construction, and a
% large shared region of near-zero pulls every correlation towards each other.
% Set data.anaMask to change it; flow_checkData derives it from ftprnt by default.
vMask = data.anaMask;
assert(isequal(size(vMask), [nrow ncol]), ...
    ['voltage mask is %s but the movie frame is %s. Run the data through ' ...
     'flow_checkData first.'], mat2str(size(vMask)), mat2str([nrow ncol]));
assert(any(vMask(:)), 'The footprint mask selects no pixel; data.ftprnt may be empty.');
fprintf('[flow] voltage correlation taken over %d anaMask pixels\n', nnz(vMask));
rVolt = get_voltCorr(maps, vRef, vMask);
rAlt  = nan(size(rVolt));
if ~isempty(vAlt), rAlt = get_voltCorr(maps, vAlt, vMask); end
score = oscScore;                  % the score is now PURELY temporal: how oscillatory the
                                   % waveform is. Where it lives is a separate, binary
                                   % question, answered by voltKill below.
voltKill = false(size(rVolt));
if c.voltPen, voltKill = abs(rVolt) >= c.voltCut; end

% ---- drop the weak rows, then order by segment ----------------------------------------
% Decided BEFORE the table is printed, so a flag beside a row is what actually happened to
% it. Both rules move across a quadrature unit, and they move in OPPOSITE directions.
scoreOK = score >= c.scoreMin;
for k = 1:nS                       % cos and sin are one quadrature UNIT: judge together
    scoreOK([k, k+nS]) = any(scoreOK([k, k+nS]));   % either half carries the pair
    % The veto propagates the other way. A quadrature pair spans every phase between
    % them, so if either half's map looks like the neuron the pair can reconstruct it --
    % rescuing the partner would let the voltage back in through the side its own map
    % does not show.
    voltKill([k, k+nS]) = any(voltKill([k, k+nS]));
end
% The veto is unconditional: dropBelowScore governs the weak rows, not this. A row that
% looks like the neuron is never regressed out, however the score knob is set.
if c.dropBelowScore, keepR = scoreOK & ~voltKill; else, keepR = ~voltKill; end
nVeto = nnz(voltKill);

fprintf('[flow] %d rows | voltage reference ''%s''\n', size(design,1), vName);
fprintf('     %-30s %8s %7s %8s\n', 'name', 'oscScore', 'r(volt)', 'score');
[~, ord] = sort(score, 'ascend');
for r = ord'
    flag = '';
    if voltKill(r)
        flag = sprintf('  <- VETOED by the voltage map (|r| >= %.2f)', c.voltCut);
    elseif ~keepR(r)
        flag = '  <- dropped';
    end
    fprintf('     %-30s %8.3f %+7.2f %8.3f%s\n', names{r}, oscScore(r), ...
        rVolt(r), score(r), flag);
end

% Recorded BEFORE the drop, because a vetoed row is not in the output to be inspected
% afterwards -- and a veto you cannot see is a veto you cannot argue with. rVoltAlt is
% the correlation against the OTHER reference: when it is much lower than rVolt, the row
% is a vessel on the dendrite rather than voltage, and the cut was too harsh.
vetoed = struct('name', {names(voltKill)}, 'rVolt', rVolt(voltKill), ...
    'rVoltAlt', rAlt(voltKill), 'oscScore', oscScore(voltKill), 'cut', c.voltCut);
if nVeto > 0
    fprintf('[flow] vetoing %d of %d rows with |r(volt)| >= %.2f\n', ...
        nVeto, numel(voltKill), c.voltCut);
end
if any(~keepR)
    nWeak = nnz(~keepR & ~voltKill);
    if nWeak > 0
        fprintf('[flow] dropping %d of %d rows below scoreMin = %.2f\n', ...
            nWeak, numel(keepR), c.scoreMin);
    end
    design = design(keepR,:);  names = names(keepR);  segOf = segOf(keepR);
    modeOf = modeOf(keepR);
    score = score(keepR);  oscScore = oscScore(keepR);  f0 = f0(keepR);
    rVolt = rVolt(keepR);  rAlt = rAlt(keepR);
    maps = maps(:,keepR);  eigImg = eigImg(:,keepR);  varExpl = varExpl(keepR);
    assert(~isempty(design), ['Every row was rejected: %d vetoed by the voltage map ' ...
        '(voltCut = %.2f), the rest below scoreMin = %.2f. Raise voltCut if the veto ' ...
        'is doing it -- but check first whether a vessel runs along the dendrite.'], ...
        nVeto, c.voltCut, c.scoreMin);
end
% ORDER BY VARIANCE EXPLAINED, largest first: the first rows in the browser are then the
% ones that actually do something to the movie, rather than whichever segment happened
% to be numbered first.
[~, os] = sort(varExpl, 'descend');
design = design(os,:);  names = names(os);  segOf = segOf(os);  modeOf = modeOf(os);
score = score(os);  oscScore = oscScore(os);  f0 = f0(os);
rVolt = rVolt(os);  rAlt = rAlt(os);
maps = maps(:,os);  eigImg = eigImg(:,os);  varExpl = varExpl(os);

% ---- what the WHOLE design removes, per pixel ----------------------------------------
% The joint fit, not the sum of the marginal ones: the rows overlap, so adding their
% individual contributions would double-count whatever they share. With Sig the design's
% own covariance and maps = cov(pixel, row), the least-squares coefficients are
% B = maps / Sig and the variance they jointly explain at a pixel is B*Sig*B', whose
% diagonal is sum(maps .* B, 2). Exact, and it needs no residual movie.
Sig  = cov(design(:,gF)');
Bj   = maps / Sig;
varFit  = reshape(sum(maps .* Bj, 2), nrow, ncol);
varMov  = reshape(vPx, nrow, ncol);
dropFrac = varFit ./ max(varMov, eps);
dropFrac(~data.anaMask) = NaN;             % outside the footprints it means nothing
onNeu = data.neuronMap > 0 & ~data.bvBW;
fprintf(['[flow] the whole design removes, inside the footprint mask: %.1f%% overall, ' ...
         '%.1f%% in the density mask,\n       %.1f%% on bvMask, %.1f%% on the DENDRITE ' ...
         'away from vessels  <- keep this one small\n'], ...
    100*mean(dropFrac(data.anaMask),'omitnan'), ...
    100*mean(dropFrac(seg.densBW & data.anaMask),'omitnan'), ...
    100*mean(dropFrac(data.bvBW & data.anaMask),'omitnan'), ...
    100*mean(dropFrac(onNeu & data.anaMask),'omitnan'));
if any(onNeu(:)) && any(data.bvBW(:)) && ...
        mean(dropFrac(onNeu),'omitnan') > 0.5*mean(dropFrac(data.bvBW & data.anaMask),'omitnan')
    warning('flow_regressors:onDendrite', ...
        ['The design removes %.1f%% on the dendrite vs %.1f%% on bvMask -- that ratio is ' ...
         'high enough that it is probably taking neural signal.'], ...
        100*mean(dropFrac(onNeu),'omitnan'), ...
        100*mean(dropFrac(data.bvBW & data.anaMask),'omitnan'));
end

reg = struct('design',design, 'names',{names}, 'segOf',segOf(:), 'modeOf',modeOf(:), ...
    'varFit',varFit, 'dropFrac',dropFrac, 'varMov',varMov, 'anaMask',data.anaMask, ...
    'maps',maps, 'eigImg',eigImg, 'varExpl',varExpl(:), ...
    'score',score(:), 'oscScore',oscScore(:), 'f0',f0(:), 'rVolt',rVolt(:), ...
    'rVoltAlt',rAlt(:), 'vetoed',vetoed, 'useSeg',useSeg, 'segPx',{segPxU}, ...
    'frames',frames, 'good',gF, 'seg',sInfo, 'srcTrace',srcNv, 'r2seg',r2, ...
    'oscCos',oscCos, 'oscSin',oscSin, 'voltRef',vRef, 'voltRefName',vName, ...
    'voltMask',vMask, 'sz',[nrow ncol]);
fprintf('[flow] design is %d x %d, segments %s\n', size(design,1), size(design,2), ...
    mat2str(useSeg(segOf)));
end

% =========================================================================
function [vRef, vName, vAlt] = voltageRef(data, cfg, fs) %#ok<INUSD>
% Which map counts as "voltage". The STA is the direct spatial signature and blood
% cannot contaminate it (the cell train is not spike-locked); the anatomical footprint is
% a proxy that also flags a vessel running ALONG the dendrite. But an STA over too few
% spikes is noise, and a noisy reference correlates with nothing -- which would switch
% the criterion off silently rather than loudly. Hence the spike-count fallback.
c = cfg.regressors;
fp = data.neuronMap ./ max(data.neuronMap(:), [], 'omitnan');
sta = [];
nT = size(data.mov_res,3);
sp = find(data.spike(:)' > 0);
sp = sp(sp > 3 & sp <= nT-3);
if numel(sp) > c.voltMinSpk
    acc = zeros(size(fp));
    for q = sp
        acc = acc + mean(-double(data.mov_res(:,:,q:q+2)),3) ...
                  - mean(-double(data.mov_res(:,:,q-3:q-1)),3);
    end
    sta = acc / numel(sp);  sta = sta / max(abs(sta(:)));
end
switch c.voltRef
    case 'auto'
        if ~isempty(sta), vRef = sta;  vName = 'sta';  vAlt = fp;
        else
            vRef = fp;  vName = 'ftprnt';  vAlt = [];
            fprintf(['[flow] voltage reference: only %d spikes (<= voltMinSpk = %d), ' ...
                     'so the STA would be noise -> ''ftprnt''\n'], numel(sp), c.voltMinSpk);
        end
    case 'sta'
        if isempty(sta)
            warning('flow_regressors:noSTA', ...
                'voltRef = ''sta'' but only %d spikes; using ''ftprnt''.', numel(sp));
            vRef = fp;  vName = 'ftprnt';  vAlt = [];
        else
            vRef = sta;  vName = 'sta';  vAlt = fp;
        end
    case 'ftprnt'
        vRef = fp;  vName = 'ftprnt';  vAlt = sta;
    otherwise
        error('cfg.regressors.voltRef must be auto | sta | ftprnt.');
end
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end
