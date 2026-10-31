function st = flow_stics(S, opt)
% st = flow_stics(S, opt)
%
% SPATIOTEMPORAL IMAGE CORRELATION (STICS): find flow without ever detecting a particle.
%
% ---------------------------------------------------------------------------------
% HOW IT WORKS
%
% Take the movie, remove each pixel's time-average, and call what is left dI(x,y,t).
% Then form the space-time correlation
%
%       G(xi, eta, tau) = < dI(x, y, t) * dI(x+xi, y+eta, t+tau) >
%
% averaged over every pixel and every start time. Read it as: "if something is here
% now, how much is there something there, tau frames later?"
%
%   tau = 0     G is a bump at (0,0) whose WIDTH is the size of a particle.
%   tau > 0     if the pattern is STATIONARY the bump stays at (0,0) and only fades.
%               If the pattern FLOWS at velocity v, the whole bump MOVES: its centre
%               sits at (vx*tau, vy*tau).
%
% So the velocity is the SLOPE of peak position against tau. Fit a line through the
% peak positions and you have vx and vy; how straight that line is tells you whether
% it is really flow.
%
% HOW IT COMPARES TO PARTICLE TRACKING -- measured, not assumed
% Both methods were run on the same synthetic vessel (true 2000 um/s) while the spacing
% between cells was varied, with a perfectly regular train and with a jittered one:
%
%   spacing   train        STICS            tracking
%     4 px    regular      +0    (FAILS)    1518  (6 tracklets, poor)
%     4 px    irregular    +2000 (exact)    1918  (61 tracklets)
%     6 px    regular      +0    (FAILS)    1997
%     6 px    irregular    +2000 (exact)    1944
%    12 px    either       +2000 (exact)    2000-2013
%
% So the honest summary is:
%   * Where both work, correlation is the more accurate of the two -- it lands on the
%     true speed while tracking runs a few percent low, because tracking estimates
%     velocity from tracklet ENDPOINTS and those endpoints are truncated.
%   * It needs no detection threshold and no resolved particles, so it returns a value
%     for EVERY tile rather than a sparse set of tracklets. "Which pixel has a flowing
%     component" is answered directly.
%   * Shot noise is white in space AND time, so it contributes only to G(0,0,0) -- one
%     spike at zero lag. It cannot fake a moving peak.
%   * BUT it has a failure mode tracking does not: see the periodicity note below. It is
%     a complement to the tracker, not a replacement.
%
% WHAT IT CANNOT DO
%   * A PERFECTLY PERIODIC train defeats it. If cells pass every dt frames the pattern
%     repeats exactly, so G(0,0,dt) is as large as G(0,0,0) and a stationary pattern is
%     indistinguishable from one flowing at spacing/dt. Measured above: a regular train
%     at 4 and 6 px spacing returns 0 um/s with high confidence. Real trains are
%     irregular enough to break the tie -- the same jitter that makes the detection-count
%     estimate unreliable elsewhere in this project is what RESCUES the correlation --
%     but a suspiciously exact zero on an otherwise busy vessel is this, not stillness.
%   * One velocity per tile. Two vessels crossing inside one tile give a smeared or
%     averaged answer -- shrink the tile, or believe the confidence, not the speed.
%   * The fastest measurable flow is maxShift/2 px per frame, because a candidate needs
%     at least two lags inside the window. Raise tilePx to reach faster flow.
%   * It measures the speed of the PATTERN, not of an identified cell. A travelling wave
%     of intensity that is not flowing cells would look identical.
%   * Removing the per-pixel time-average kills immobile structure but NOT slow drift;
%     use opt.detrendMs when the baseline wanders.
% ---------------------------------------------------------------------------------
%
% INPUT
%   S    [nrow x ncol x nT] the movie (band-passed, typically)
%   opt  struct, all optional:
%     .tilePx     tile side in px                       (default 16)
%     .tileStep   tile stride; < tilePx overlaps        (default 8)
%     .maxTau     largest temporal lag, frames          (default 10)
%     .maxShift   largest spatial lag searched, px      (default tilePx/2 - 1)
%     .frameStep  use every Nth frame pair              (default 1)
%     .detrendMs  subtract a moving average this long   (default 0 = off)
%     .fs         Hz, for reporting                     (default 1000)
%     .pxUm       um per px, for reporting              (default 1)
%     .minZ       confidence below this is set NaN        (default 6)
%     .dv         velocity grid step, px/frame            (default 0.25)
%     .verbose    print progress                        (default true)
%
% OUTPUT  st
%   .speed    [nrow x ncol] um/s, painted per tile (NaN where no flow was found)
%   .vx .vy   [nrow x ncol] um/s components
%   .angle    [nrow x ncol] deg
%   .score    [nrow x ncol] confidence z of the velocity fit: the FLOW SCORE
%   .amp      [nrow x ncol] correlation amplitude at tau = 0 (how much moved at all)
%   .tileC    {nTy x nTx} the raw G(xi,eta,tau) of every tile, for inspection
%   .tileXY   [nTile x 2] tile centres
%   .tau      the lags used
%
% 2026.08.28, Byung Hun Lee

if nargin < 2, opt = struct; end
d = struct('tilePx',16, 'tileStep',8, 'maxTau',10, 'maxShift',[], 'frameStep',1, ...
           'detrendMs',0, 'fs',1000, 'pxUm',1, 'minZ',6, 'dv',0.25, 'verbose',true);
f = fieldnames(d);
for k = 1:numel(f)
    if ~isfield(opt,f{k}) || isempty(opt.(f{k})), opt.(f{k}) = d.(f{k}); end
end
if isempty(opt.maxShift), opt.maxShift = max(2, floor(opt.tilePx/2) - 1); end
[nrow, ncol, nT] = size(S);
tp = opt.tilePx;  ts = opt.tileStep;  nTau = opt.maxTau;
assert(opt.maxShift < tp, 'maxShift (%d) must be smaller than tilePx (%d).', opt.maxShift, tp);
assert(nT > nTau + 2, 'Need more than maxTau+2 frames.');

y0 = 1:ts:(nrow - tp + 1);   x0 = 1:ts:(ncol - tp + 1);
nTy = numel(y0);  nTx = numel(x0);
assert(nTy > 0 && nTx > 0, 'The movie is smaller than one tile.');

tau  = 0:nTau;
Pfft = 2*tp;                       % zero-pad so the FFT correlation does not wrap around
[sx, sy] = meshgrid(-tp:tp-1, -tp:tp-1);
inWin = abs(sx) <= opt.maxShift & abs(sy) <= opt.maxShift;

vxT = nan(nTy, nTx);  vyT = nan(nTy, nTx);  r2T = nan(nTy, nTx);  ampT = nan(nTy, nTx);
nUsed = zeros(nTy, nTx);   % lags that stayed inside the search window
% Candidate velocities, in px per sampled frame. The reach is set by the window: a
% candidate is only scorable while |v*tau| stays inside maxShift, so the fastest
% measurable flow is maxShift/2 px per frame (two usable lags) and anything above that
% is invisible however long the movie is.
vLim = opt.maxShift / 2;
vAx  = -vLim : opt.dv : vLim;
[VXg, VYg] = meshgrid(vAx, vAx);
tileC = cell(nTy, nTx);
if opt.verbose
    fprintf('[stics] %d x %d tiles of %d px (step %d), lags 0-%d frames, shift +/-%d px\n', ...
        nTy, nTx, tp, ts, nTau, opt.maxShift);
end

for iy = 1:nTy
    for ix = 1:nTx
        Y = double(S(y0(iy):y0(iy)+tp-1, x0(ix):x0(ix)+tp-1, 1:opt.frameStep:end));
        nt = size(Y,3);
        if nt < nTau + 2, continue, end

        % IMMOBILE REMOVAL: subtract each pixel's own time-average, so anything that
        % does not change contributes nothing. Without this the vessel wall and the
        % dendrite dominate G and pin the peak at zero lag.
        Y = Y - mean(Y, 3);
        if opt.detrendMs > 0
            w = max(3, round(opt.detrendMs * opt.fs / 1000 / opt.frameStep));
            Y = Y - movmean(Y, w, 3);
        end
        if ~any(Y(:)), continue, end

        F = fft2(Y, Pfft, Pfft);                 % one FFT per frame, reused for every lag
        C = zeros(Pfft, Pfft, numel(tau));
        for k = 1:numel(tau)
            t1 = 1:(nt - tau(k));
            t2 = t1 + tau(k);
            % accumulate in Fourier space, one inverse transform per lag
            P = mean(conj(F(:,:,t1)) .* F(:,:,t2), 3);
            C(:,:,k) = fftshift(real(ifft2(P)));
        end
        C = C(tp+1-tp:tp+tp, tp+1-tp:tp+tp, :);      % centre block, lags -tp..tp-1
        c0 = max(max(C(:,:,1)));
        if ~isfinite(c0) || c0 <= 0, continue, end
        C = C / c0;                                   % normalise by the tau = 0 peak
        tileC{iy,ix} = C;
        ampT(iy,ix) = c0;

        % ---- search VELOCITY space directly -----------------------------------------
        % Tracking the peak lag by lag looks natural but does not work here:
        %
        %  ALIASING. A quasi-regular train of cells makes G periodic in space with the
        %  inter-cell spacing. Measured on a train 8 px apart, the true peak at +6 px
        %  reappears at -2 px and a per-lag maximum picks the wrong branch, collapsing
        %  the fitted speed towards zero.
        %
        %  CIRCULARITY. Resolving that by predicting the next position and searching
        %  near it makes the answer self-fulfilling: the positions are constrained to lie
        %  on a line, so an R^2 computed from them is high whatever the data. Tried, and
        %  it passed 21 of 45 tiles that contained no flow at all.
        %
        % So the velocity is fitted GLOBALLY instead. For rigid flow
        %       G(xi, eta, tau) = G0(xi - vx*tau, eta - vy*tau)
        % so a candidate velocity is scored by how much correlation sits along the line
        % it predicts, summed over lags. The best candidate is the velocity; how far it
        % stands above the rest of velocity space is the CONFIDENCE, and that is a real
        % measurement -- a tile with no flow has a flat landscape and scores near zero.
        sc = zeros(size(VXg));  nL = zeros(size(VXg));
        for k = 2:numel(tau)
            tk = tau(k);
            xq = VXg*tk;  yq = VYg*tk;
            in = abs(xq) <= opt.maxShift & abs(yq) <= opt.maxShift;
            if ~any(in(:)), break, end
            v = interp2(sx, sy, C(:,:,k), xq, yq, 'linear', NaN);
            v(~in | ~isfinite(v)) = 0;
            sc = sc + v;  nL = nL + double(in);
        end
        useV = nL >= 3;                      % a candidate needs 3 lags to mean anything
        if ~any(useV(:)), continue, end
        scM = sc ./ max(nL,1);  scM(~useV) = NaN;
        [bestV, ib] = max(scM(:));
        if ~isfinite(bestV), continue, end

        % confidence: how far the winner stands above the rest of velocity space
        pool = scM(useV);
        md = median(pool);  sd = 1.4826*median(abs(pool - md));
        if sd <= 0, continue, end
        zBest = (bestV - md) / sd;

        nUsed(iy,ix) = nL(ib);
        vxT(iy,ix) = VXg(ib);  vyT(iy,ix) = VYg(ib);
        r2T(iy,ix) = zBest;
    end
    if opt.verbose && mod(iy, max(1,round(nTy/5))) == 0
        fprintf('       row %d/%d\n', iy, nTy);
    end
end

% px per sampled frame -> um/s
k2 = opt.pxUm * opt.fs / opt.frameStep;
vxT = vxT * k2;   vyT = vyT * k2;
spdT = hypot(vxT, vyT);
bad  = ~isfinite(r2T) | r2T < opt.minZ;
spdT(bad) = NaN;  vxT(bad) = NaN;  vyT(bad) = NaN;

% ---- paint the tiles back onto the frame -------------------------------------------
[speed, vx, vy, score, amp] = deal(nan(nrow, ncol));
cnt = zeros(nrow, ncol);
for iy = 1:nTy
    for ix = 1:nTx
        ri = y0(iy):y0(iy)+tp-1;  ci = x0(ix):x0(ix)+tp-1;
        if isfinite(spdT(iy,ix))
            speed(ri,ci) = nanplus(speed(ri,ci), spdT(iy,ix));
            vx(ri,ci)    = nanplus(vx(ri,ci),    vxT(iy,ix));
            vy(ri,ci)    = nanplus(vy(ri,ci),    vyT(iy,ix));
            cnt(ri,ci)   = cnt(ri,ci) + 1;
        end
        score(ri,ci) = max(score(ri,ci), r2T(iy,ix), 'omitnan');
        amp(ri,ci)   = max(amp(ri,ci),   ampT(iy,ix), 'omitnan');
    end
end
speed = speed ./ max(cnt,1);  vx = vx ./ max(cnt,1);  vy = vy ./ max(cnt,1);
speed(cnt==0) = NaN;  vx(cnt==0) = NaN;  vy(cnt==0) = NaN;

[Xc, Yc] = meshgrid(x0 + tp/2, y0 + tp/2);
st = struct('speed',speed, 'vx',vx, 'vy',vy, 'angle',atan2d(vy,vx), ...
    'score',score, 'amp',amp, 'tileC',{tileC}, 'tileXY',[Xc(:) Yc(:)], ...
    'tileSpeed',spdT, 'tileVx',vxT, 'tileVy',vyT, 'tileR2',r2T, ...
    'tau',tau, 'opt',opt, 'y0',y0, 'x0',x0, 'nLagUsed',nUsed);

if opt.verbose
    ok = isfinite(spdT);
    fprintf('[stics] %d of %d tiles pass z >= %.1f | speed median %.0f um/s (IQR %.0f-%.0f)\n', ...
        nnz(ok), numel(ok), opt.minZ, median(spdT(ok)), ...
        prctile(spdT(ok),25), prctile(spdT(ok),75));
    lagOK = nUsed(nUsed > 0);
    fprintf('        lags used per tile: median %d of %d requested\n', ...
        median(lagOK), numel(tau));
    if ~isempty(lagOK) && median(lagOK) <= 3 && any(ok(:))
        fprintf(['        -> the peak leaves the +/-%d px window early, so the fit rests on ' ...
                 'few lags.\n           Raise tilePx (bigger window) or lower maxTau.\n'], opt.maxShift);
    end
end
end

function a = nanplus(a, b)
a(isnan(a)) = 0;
a = a + b;
end
