function ap = flow_apply(data, bp, seg, reg, cfg)
% ap = flow_apply(data, bp, seg, reg, cfg)
%
% Take the design that flow_regressors built on the tracked window and apply it to the
% WHOLE recording: extend it in time, fit it, and measure what it removes.
%
% =====================================================================================
% WHAT CARRIES OVER AND WHAT DOES NOT
%
% Two different objects came out of the earlier stages and they behave differently:
%
%   SPATIAL, reusable.  seg.px (which pixels are a vessel segment) and reg.seg(k).U (the
%   spatial modes inside each). These say WHERE the vessel is. The movie is motion-
%   corrected and mean-subtracted, so the vessel is in the same place at the last frame
%   as at the first; re-running the detection on the rest would re-estimate the same mask
%   from a different noise realisation.
%
%   TEMPORAL, must be recomputed.  The oscillators are A(t)cos(phi(t)), defined only on
%   the frames they were built from. They are rebuilt here by projecting the SAME modes
%   onto the full record.
%
%   THE PASSBAND can go stale. fPeak was measured on the tracked window; if flow speeds
%   up or slows down the peak moves. cfg.apply.bandMode = 'refit' re-finds it, and the
%   per-block diagnostic below shows whether it drifted enough to matter.
%
% refitPCA goes further and re-derives the modes themselves on the full record. The masks
% are geometry and do not change; the modes are a STATISTICAL estimate from a short
% window, and are the one reused quantity worth redoing on more data. When it is on, the
% row selection is redone too, because flow_regressors is what runs.
%
% THE TIME-VARYING-COEFFICIENT TEST. A single beta per pixel assumes the vessel's
% contribution is constant for the whole recording. cfg.apply.split refits it in
% overlapping windows and asks whether that is better -- OUT OF SAMPLE, because a
% windowed fit has ~nWin times the free parameters and therefore wins in-sample always.
% Whole blocks are held out, not single frames, since a lone held-out frame is
% predictable from its neighbours and would hand the flexible model a free pass.
% =====================================================================================
%
% OUTPUT  ap
%   .design .names .segOf     the full-record design
%   .C        [nPix x (1+nReg)] fit coefficients, column 1 = DC. The corrected movie is
%             mov - C*[1; centred design], and any trace follows in closed form as
%             T_corr = T_raw + (F'*C)*I -- no corrected movie is ever built.
%   .varIn .varFit .dropFrac  per-pixel variance in, removed, and the ratio
%   .fPeak .band .drift       per segment, and how far the peak moved across the record
%   .split    the held-out comparison, when cfg.apply.split is on. .varFitG and .varFitS
%             are var(fitted) for the global and the windowed fit, both measured on the
%             SAME held-out frames -- the only footing on which they are comparable.
%
% 2026.08.28, Byung Hun Lee

c  = cfg.apply;
fs = cfg.common.fs;
[nrow, ncol, nT] = size(data.mov_res);
npx = nrow*ncol;
gAll = ~data.badfr(:)';  nG = nnz(gAll);
fFloor = max(pick(c.fMin, cfg.regressors.fMin), 1);
poolGuard = flow_noPool();  %#ok<NASGU>

fprintf('[flow] apply: extending %d regressors from %d frames to all %d (%.1f s, %d good)\n', ...
    size(reg.design,1), numel(reg.frames), nT, nT/fs, nG);

% ---- rebuild the design over every frame ---------------------------------------------
if c.refitPCA
    % Re-derive everything on the full record. flow_regressors already does exactly this
    % given a different frame range, so it is called again rather than duplicated -- the
    % two paths can then never drift apart.
    fprintf('[flow] apply: re-deriving the modes on the full record (refitPCA)\n');
    trkFull = struct('frames', 1:nT);
    regF = flow_regressors(bp, trkFull, seg, data, cfg);
    design = regF.design;  names = regF.names;  segOf = regF.segOf;
    sInfo  = regF.seg;     useSeg = regF.useSeg;
else
    [design, names, segOf, sInfo, useSeg] = extendDesign(bp, reg, cfg, fs, nT, gAll, fFloor);
    fprintf('[flow] apply: reused the (step 4) modes, %d rows rebuilt over %d frames\n', ...
        size(design,1), nT);
end
nReg = size(design,1);

% ---- did the band drift? --------------------------------------------------------------
% A block cannot be longer than the record. Without the clamp, floor(nT/blk) is 0 for a
% short recording, max(1,0) forces one block anyway, and that block still indexes blk
% frames past the end.
blk = max(round(c.blockS*fs), 64);
if blk > nT
    fprintf(['[flow] apply: blockS = %g s is %d frames but the record is %d; using the ' ...
             'whole record as one block (no drift estimate is possible)\n'], ...
        c.blockS, blk, nT);
    blk = nT;
end
nBlk = max(1, floor(nT/blk));
fPk = nan(numel(sInfo), nBlk);
for k = 1:numel(sInfo)
    if ~isfield(sInfo, 'srcFull') || isempty(sInfo(k).srcFull), continue, end
    y = sInfo(k).srcFull;
    if numel(y) ~= nT, continue, end        % built on a different frame range
    for b = 1:nBlk
        ib = (b-1)*blk + (1:blk);
        yy = y(ib) - mean(y(ib));  yy(~gAll(ib)) = 0;
        w  = 0.5*(1 - cos(2*pi*(0:blk-1)/max(blk-1,1)));
        P  = abs(fft(yy.*w)).^2;  fv = (0:blk-1)*fs/blk;
        in = fv >= sInfo(k).srch(1) & fv <= sInfo(k).srch(2) & fv < fs/2;
        if any(in), P(~in) = -inf; [~,ip] = max(P); fPk(k,b) = fv(ip); end
    end
end
drift = (max(fPk,[],2,'omitnan') - min(fPk,[],2,'omitnan')) ./ ...
        max(arrayfun(@(x) x.fPeak, sInfo(:)), eps);
fprintf('[flow] band drift across %d blocks of %.2f s: median %.2f, max %.2f\n', ...
    nBlk, blk/fs, median(drift,'omitnan'), max(drift,[],'omitnan'));
if any(drift > cfg.regressors.bandFrac)
    warning('flow_apply:drift', ...
        ['%d segment(s) move their spectral peak by more than bandFrac = %.2f across the ' ...
         'record: one band cannot cover the whole movie there.'], ...
        nnz(drift > cfg.regressors.bandFrac), cfg.regressors.bandFrac);
end

% ---- fit, without ever copying the movie ----------------------------------------------
% Same normal equations SeeResiduals solves, accumulated in time chunks. The per-pixel
% mean and variance of the input come out of the same pass.
switch lower(c.fitOn)
    case 'mov_res', src = data.mov_res;   srcNm = 'mov_res (voltage still in it)';
    case 'mov_nv',  src = bp.mov_nv;      srcNm = 'mov_nv (voltage regressed out)';
    otherwise, error('cfg.apply.fitOn must be ''mov_res'' or ''mov_nv''.');
end
fprintf('[flow] apply: fitting to %s on %d good frames\n', srcNm, nG);
I  = [ones(1,nG); design(:,gAll) - mean(design(:,gAll),2)];
A  = I*I';
gi = find(gAll);
M  = zeros(npx, nReg+1);  s1 = zeros(npx,1);  s2 = zeros(npx,1);
for c0 = 1:c.fitChunk:nG
    jj = c0:min(c0+c.fitChunk-1, nG);
    Y  = double(reshape(src(:,:,gi(jj)), npx, numel(jj)));
    M  = M + Y*I(:,jj)';
    s1 = s1 + sum(Y,2);  s2 = s2 + sum(Y.^2,2);
end
clear Y
C = M / A;
varIn  = reshape((s2 - s1.^2/nG)/max(nG-1,1), nrow, ncol);
Sg     = cov(design(:,gAll)');
varFit = reshape(sum((C(:,2:end)*Sg).*C(:,2:end), 2), nrow, ncol);
dropFrac = varFit ./ max(varIn, eps);
dropFrac(~data.anaMask) = NaN;      % every fraction below is inside the footprint mask
onNeu = data.neuronMap > 0 & ~data.bvBW;
mm = @(m) 100*mean(dropFrac(m & data.anaMask), 'omitnan');
fprintf(['[flow] removed inside the footprint mask (%d px): %.1f%% overall, %.1f%% in the ' ...
         'density mask,\n       %.1f%% on bvMask, %.1f%% on the DENDRITE away from ' ...
         'vessels  <- keep this one small\n'], nnz(data.anaMask), ...
    100*mean(dropFrac(data.anaMask),'omitnan'), mm(seg.densBW), mm(data.bvBW), mm(onNeu));

ap = struct('design',design, 'names',{names}, 'segOf',segOf(:), 'C',C, ...
    'varIn',varIn, 'varFit',varFit, 'dropFrac',dropFrac, 'anaMask',data.anaMask, 'good',gAll, ...
    'seg',sInfo, 'useSeg',useSeg, 'fPeakBlk',fPk, 'drift',drift, ...
    'fitOn',c.fitOn, 'sz',[nrow ncol], 'split',[]);

% ---- is a time-varying coefficient better? --------------------------------------------
if c.split
    ap.split = splitTest(src, design, gAll, c, nrow, ncol, npx, ...
                         data, seg, varIn);
end
end

% =========================================================================
function [design, names, segOf, sInfo, useSeg] = extendDesign(bp, reg, cfg, fs, nT, gAll, fFloor)
% Project the STORED modes onto the full record. Complex modes must be projected onto the
% ANALYTIC signal, not the raw traces: they were derived from hilbert(bandpass(Y)), so
% projecting the real traces would mix the two quadratures and lose the phase the mode
% encodes.
c = cfg.regressors;
Sv = reshape(bp.S, [], size(bp.S,3));
nS = numel(reg.useSeg);
oscCos = zeros(nS,nT); oscSin = oscCos;
sInfo = reg.seg;
for k = 1:nS
    px = reg.segPx{k};
    Y  = -double(Sv(px,:));  Y = Y - mean(Y,2);
    U  = sInfo(k).U;
    if ~isempty(U)
        if c.cpca
            Yb = get_bandpass(Y, fs, sInfo(k).srch);  Yb(~isfinite(Yb)) = 0;
            Z  = hilbert(Yb.').';  Z = Z - mean(Z,2);
            pc = U' * Z;  y = real(pc(1,:));
        else
            pc = U' * Y;  y = pc(1,:);
        end
    else
        pc = [];  y = mean(Y,1);
    end
    sInfo(k).pcAll   = pc;
    sInfo(k).srcFull = y;

    if strcmpi(cfg.apply.bandMode,'reuse')
        fPeak = max(sInfo(k).fPeak, fFloor);
    else
        yg = y;  yg(~gAll) = NaN;
        [Py,Fy] = nanPSD(yg, fs, min(1024, 2^floor(log2(nT/2))));
        Py = Py(:); Fy = Fy(:);
        in = Fy >= sInfo(k).srch(1) & Fy <= sInfo(k).srch(2);
        if ~any(in & isfinite(Py)), in = Fy > fFloor & Fy < fs/2-1; end
        Pm = Py; Pm(~in | ~isfinite(Pm)) = -inf;
        [~,ip] = max(Pm); fPeak = Fy(ip);
    end
    bnd = fPeak*[1-c.bandFrac, 1+c.bandFrac];
    bnd = [max(bnd(1),fFloor), min(bnd(2), fs/2-1)];    % floor AFTER the widening
    if ~all(isfinite(bnd)) || bnd(2) <= bnd(1)*1.05, bnd = sInfo(k).srch; end
    sInfo(k).fPeak = fPeak;  sInfo(k).band = bnd;

    yb = get_bandpass(y, fs, bnd); yb(~isfinite(yb)) = 0;
    z = hilbert(yb);
    oscCos(k,:) = abs(z).*cos(angle(z));
    oscSin(k,:) = abs(z).*sin(angle(z));
end
zrow = @(V)(V-mean(V,2))./max(std(V,0,2),eps);
% same canonical order as step 4, so the rows mean the same thing
base = cell(1,nS);
for k = 1:nS, base{k} = strtok(reg.names{find(reg.segOf==k,1)}, '|'); end
design = [zrow(oscCos); zrow(oscSin)];
names  = [cellfun(@(v)[strtrim(v) ' | cos'], base,'uni',0), ...
          cellfun(@(v)[strtrim(v) ' | sin'], base,'uni',0)];
segOf  = [1:nS, 1:nS];
for k = 1:nS
    P = sInfo(k).pcAll;
    if isempty(P) || size(P,1) < 2, continue, end
    ex = P(2:end,:);  q = 2:size(P,1);
    if c.cpca
        design = [design; zrow(real(ex)); zrow(imag(ex))];                 %#ok<AGROW>
        names  = [names, arrayfun(@(j)sprintf('%s | extra%d re',strtrim(base{k}),j),q,'uni',0), ...
                         arrayfun(@(j)sprintf('%s | extra%d im',strtrim(base{k}),j),q,'uni',0)]; %#ok<AGROW>
        segOf  = [segOf, repmat(k,1,2*size(ex,1))];                        %#ok<AGROW>
    else
        design = [design; zrow(ex)];                                       %#ok<AGROW>
        names  = [names, arrayfun(@(j)sprintf('%s | extra%d',strtrim(base{k}),j),q,'uni',0)]; %#ok<AGROW>
        segOf  = [segOf, repmat(k,1,size(ex,1))];                          %#ok<AGROW>
    end
end
useSeg = reg.useSeg;
end

% =========================================================================
function sp = splitTest(src, design, gAll, c, nrow, ncol, npx, data, seg, varIn) %#ok<INUSD>
% Global vs windowed coefficients, scored on held-out blocks.
nT = size(design,2);  nReg = size(design,1);
win = min(c.splitWin, nT);              % likewise: a window cannot exceed the record
ovl = min(c.splitOvl, max(win-1,0));
if win < c.splitWin
    fprintf(['[flow] split: splitWin %d > %d frames available; using %d with %d ' ...
             'overlap\n'], c.splitWin, nT, win, ovl);
end
stepW = max(1, win - ovl);
wS = 1:stepW:max(1, nT - win + 1);
if wS(end) + win - 1 < nT, wS(end+1) = nT - win + 1; end
wS = unique(max(wS,1));  nWin = numel(wS);  wE = min(wS + win - 1, nT);
W = zeros(nWin, nT);
for k = 1:nWin
    L = wE(k)-wS(k)+1;  tp = ones(1,L);  nr = min(ovl, floor(L/2));
    if nr > 0
        ramp = (1:nr)/(nr+1);
        if k > 1,    tp(1:nr) = ramp; end
        if k < nWin, tp(end-nr+1:end) = fliplr(ramp); end
    end
    W(k, wS(k):wE(k)) = tp;
end
W = W ./ sum(W,1);
chunkId = ceil((1:nT)/c.cvChunk);
testM  = mod(chunkId, max(2,round(1/max(c.cvFrac,eps)))) == 0;
trainF = gAll & ~testM;  testF = gAll & testM;
fprintf('[flow] split test: %d windows of %d (overlap %d), train %d, held out %d\n', ...
    nWin, win, ovl, nnz(trainF), nnz(testF));
if nWin < 2
    fprintf(['       only one window fits, so this cannot test anything -- lower ' ...
             'cfg.apply.splitWin or skip the split test\n']);
end

I = [ones(1,nT); design - mean(design(:,trainF),2)];
Cg = fitOn(src, I, find(trainF), c.fitChunk, npx, nReg);
Cw = cell(1,nWin);  nFall = 0;
for k = 1:nWin
    idx = find(trainF & (1:nT) >= wS(k) & (1:nT) <= wE(k));
    if numel(idx) < c.splitMinR*(nReg+1), Cw{k} = Cg; nFall = nFall+1; continue, end
    Cw{k} = fitOn(src, I, idx, c.fitChunk, npx, nReg);
end
if nFall > 0
    fprintf('       %d of %d windows fell back to the global coefficients\n', nFall, nWin);
end
te = find(testF);  nTe = numel(te);
s = zeros(npx,3); q = zeros(npx,3);
for c0 = 1:c.fitChunk:nTe
    jj = c0:min(c0+c.fitChunk-1, nTe);  ci = te(jj);
    Y  = single(reshape(src(:,:,ci), npx, numel(ci)));
    Ic = single(I(:,ci));
    Rg = Y - single(Cg)*Ic;
    Rs = Y;
    for k = 1:nWin
        act = find(W(k,ci) > 0);
        if isempty(act), continue, end
        Rs(:,act) = Rs(:,act) - single(Cw{k})*(Ic(:,act).*single(W(k,ci(act))));
    end
    s(:,1)=s(:,1)+sum(double(Y),2);  q(:,1)=q(:,1)+sum(double(Y).^2,2);
    s(:,2)=s(:,2)+sum(double(Rg),2); q(:,2)=q(:,2)+sum(double(Rg).^2,2);
    s(:,3)=s(:,3)+sum(double(Rs),2); q(:,3)=q(:,3)+sum(double(Rs).^2,2);
end
v = (q - s.^2/nTe) / max(nTe-1,1);
% var(raw) - var(residual) IS the variance the fit removed, so both fitted variances come
% out of the same held-out pass with no extra work and no fitted movie.
varFitG = reshape(max(v(:,1)-v(:,2), 0), nrow, ncol);
varFitS = reshape(max(v(:,1)-v(:,3), 0), nrow, ncol);
varRawT = reshape(v(:,1), nrow, ncol);
varFitG(~data.anaMask) = NaN;  varFitS(~data.anaMask) = NaN;
dropG = reshape(1 - v(:,2)./max(v(:,1),eps), nrow, ncol);
dropS = reshape(1 - v(:,3)./max(v(:,1),eps), nrow, ncol);
dropG(~data.anaMask) = NaN;  dropS(~data.anaMask) = NaN;
onNeu = data.neuronMap > 0 & ~data.bvBW;
gV = 100*mean(dropS(data.bvBW & data.anaMask) - dropG(data.bvBW & data.anaMask), 'omitnan');
gD = 100*mean(dropS(onNeu & data.anaMask)     - dropG(onNeu & data.anaMask), 'omitnan');
fprintf('       HELD-OUT removed (%%):     global   split   split-global\n');
mk = {data.bvBW,'bvMask'; seg.densBW,'density'; onNeu,'dendrite'};
for k = 1:3
    m = mk{k,1} & data.anaMask;
    fprintf('         %-10s %14.2f %7.2f %+12.2f\n', mk{k,2}, ...
        100*mean(dropG(m),'omitnan'), 100*mean(dropS(m),'omitnan'), ...
        100*mean(dropS(m)-dropG(m),'omitnan'));
end
if gV > 0.5 && gD < 0.5*gV
    verdict = sprintf(['SPLIT WINS: +%.2f%% more on the vessels out of sample while the ' ...
        'dendrite moves +%.2f%%. The coefficient really is time-varying.'], gV, gD);
elseif gV <= 0
    verdict = sprintf(['GLOBAL WINS: the windowed fit removes no more on the vessels out ' ...
        'of sample (%+.2f%%), so the extra parameters only fit noise.'], gV);
else
    verdict = sprintf(['AMBIGUOUS: +%.2f%% on the vessels but +%.2f%% on the dendrite -- ' ...
        'the windowed fit is buying vessel variance partly with neural signal.'], gV, gD);
end
fprintf('       -> %s\n', verdict);
bvIdx = find(data.bvBW);
betaW = zeros(nWin, nReg);
for k = 1:nWin, betaW(k,:) = mean(Cw{k}(bvIdx,2:end), 1); end
sp = struct('dropG',dropG, 'dropS',dropS, 'varFitG',varFitG, 'varFitS',varFitS, ...
    'varRaw',varRawT, 'betaW',betaW, ...
    'betaG',mean(Cg(bvIdx,2:end),1), 'wStart',wS, 'wEnd',wE, 'W',W, ...
    'testMask',testM, 'gainVessel',gV, 'gainDendrite',gD, 'verdict',verdict);
end

function C = fitOn(src, I, idx, chunk, npx, nReg)
A = zeros(nReg+1);  M = zeros(npx, nReg+1);
for c0 = 1:chunk:numel(idx)
    jj = c0:min(c0+chunk-1, numel(idx));  ci = idx(jj);
    Y = double(reshape(src(:,:,ci), npx, numel(ci)));
    M = M + Y*I(:,ci)';
    A = A + I(:,ci)*I(:,ci)';
end
if rcond(A) < 1e-12, C = zeros(npx, nReg+1); else, C = M / A; end
end

function v = pick(a, b)
if isempty(a), v = b; else, v = a; end
end
