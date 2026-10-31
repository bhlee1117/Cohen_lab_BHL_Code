function cp = flow_compare(data, ap, reg, cfg)
% cp = flow_compare(data, ap, reg, cfg)
%
% What the correction does to the TRACES -- the quantities that actually get reported.
%
% NO CORRECTED MOVIE IS BUILT. The correction is a linear projection, so for any
% footprint F the corrected trace follows in closed form:
%
%       T_corr = -F'*(mov - C*I) = T_raw + (F'*C)*I
%
% with C = ap.C and I = [1; centred design]. F'*C is nROI x (1+nReg), so the whole thing
% is one small matrix product instead of a second copy of an 8 GB movie.
%
% THREE CHECKS, in increasing order of how much they matter:
%   1. WHERE the removed variance sits. It must concentrate on ROIs that overlap bvMask.
%      A correction that removes as much off the vessels as on them is not vascular.
%   2. WHICH FREQUENCIES it comes out of. The regressors are narrowband, so the before/
%      after PSDs should differ inside the segment bands and lie on top of each other
%      elsewhere. A broadband drop means the design is soaking up general variance.
%   3. THE SPIKE-TRIGGERED AVERAGE. Spike frames were excluded from the fit, so the
%      REMOVED component's own STA is a genuine out-of-sample test and should be flat.
%      If it has a deflection at t = 0 the regressor has learned the AP.
%
% .segFrac records how much of each ROI lies on a regressor segment, so the removed
% variance can be read against where the design was actually derived from.
%
% NORMALISATION. Traces are divided by the bAP peak measured on their OWN spike-triggered
% average, so 1.0 = one spike and ROIs of different depth are comparable. Measuring the
% bAP from the trace itself is what makes this independent of how the trace was extracted
% -- these project onto the UNNORMALISED footprints, so any externally fitted scale
% factor would be off by an unknown gain.
%
% 2026.08.28, Byung Hun Lee

c  = cfg.compare;
fs = cfg.common.fs;
Wsm = cfg.common.smoothMs;
[nrow, ncol, nT] = size(data.mov_res);
npx = nrow*ncol;
gd  = ~data.badfr(:)';
assert(~isempty(data.ftprnt), 'flow_compare needs data.ftprnt.');

% ---- traces, before and after, in closed form ----------------------------------------
% Footprints are used UNNORMALISED, so the only difference between the traces being
% compared is the movie behind them and not the weighting. Scale is not a problem: the
% normalisation below divides by each trace's own measured bAP.
F = reshape(double(data.ftprnt), npx, []);
nR = size(F,2);
I  = [ones(1,nT); ap.design - mean(ap.design(:,ap.good), 2)];
Traw = -(F' * reshape(data.mov_res, npx, nT));      % spike-positive
Tbld = -(F' * ap.C) * I;                            % the fitted vascular component
Tcor = Traw - Tbld;

dfb  = @(T) T - movmedian(T, c.baseWin, 2, 'omitnan');
dRaw = dfb(Traw);  dCor = dfb(Tcor);  dBld = dfb(Tbld);
fRaw = dRaw - movmean(dRaw, Wsm, 2);   fCor = dCor - movmean(dCor, Wsm, 2);
sRaw = movmean(dRaw, Wsm, 2);          sCor = movmean(dCor, Wsm, 2);

% ---- 1. where does the removed variance sit? -----------------------------------------
ftBin = reshape(double(data.ftprnt) > 0, npx, nR);
npxR  = sum(ftBin, 1)';
bvFrac = (ftBin' * double(data.bvBW(:))) ./ max(npxR,1);
onVes  = bvFrac >= c.bvFrac;

% Which ROIs actually sit on a segment the design was built from. This is the honest
% expectation: those are the ROIs the correction should touch, and any OTHER ROI losing
% variance is the correction reaching somewhere it was not derived from.
segMask = false(nrow, ncol);
if isfield(reg,'segPx')
    for k = 1:numel(reg.segPx), segMask(reg.segPx{k}) = true; end
end
segFrac = (ftBin' * double(segMask(:))) ./ max(npxR,1);
onSeg   = segFrac > 0;
fprintf('[flow] compare: %d of %d ROIs overlap a regressor segment\n', nnz(onSeg), nR);
remAll  = 100*(1 - var(dCor(:,gd),0,2)./max(var(dRaw(:,gd),0,2),eps));
remFast = 100*(1 - var(fCor(:,gd),0,2)./max(var(fRaw(:,gd),0,2),eps));
remSub  = 100*(1 - var(sCor(:,gd),0,2)./max(var(sRaw(:,gd),0,2),eps));
fprintf('[flow] compare: %d ROIs, %d with >= %.0f%% of their pixels on bvMask\n', ...
    nR, nnz(onVes), 100*c.bvFrac);
fprintf('       variance removed (median over ROIs):  on vessel   off vessel   ratio\n');
pr = @(nm,v) fprintf('         %-20s %9.2f %12.2f %7.1f\n', nm, ...
    median(v(onVes)), median(v(~onVes)), median(v(onVes))/max(median(v(~onVes)),eps));
if any(onVes) && any(~onVes)
    pr('full band (%)', remAll);  pr(sprintf('fast <%d ms (%%)',Wsm), remFast);
    pr(sprintf('%d ms subthr (%%)',Wsm), remSub);
    if median(remAll(~onVes)) > 0.5*median(remAll(onVes))
        warning('flow_compare:offVessel', ...
            ['The correction removes nearly as much off the vessels as on them ' ...
             '(%.2f%% vs %.2f%%) -- it is not behaving like a vascular regressor.'], ...
            median(remAll(~onVes)), median(remAll(onVes)));
    end
else
    fprintf('         (every ROI falls on the same side of bvFrac; ratio not computable)\n');
end

% ---- 3. the spike-triggered average --------------------------------------------------
w   = c.staWin(1):c.staWin(2);
sIx = find(data.spike(:)' > 0);
sIx = sIx(sIx + w(1) >= 1 & sIx + w(end) <= nT);
sta = struct('w',w,'raw',[],'cor',[],'bld',[],'nSpk',numel(sIx), ...
             'ampRaw',[],'ampCor',[],'ampBld',[],'bapMeas',[]);
if numel(sIx) >= 5
    idxV = reshape(sIx(:)' + w(:), [], 1);  nW = numel(w);
    bB = w >= c.staBase(1) & w <= c.staBase(2);
    bP = w >= c.staPeak(1) & w <= c.staPeak(2);
    one = @(D,r) mean(reshape(D(r, idxV), nW, []), 2)';
    [sr, sc_, sb] = deal(nan(nR, nW));
    for r = 1:nR
        a = one(dRaw,r); sr(r,:) = a - mean(a(bB));
        a = one(dCor,r); sc_(r,:) = a - mean(a(bB));
        a = one(dBld,r); sb(r,:) = a - mean(a(bB));
    end
    sta.raw = sr; sta.cor = sc_; sta.bld = sb;
    sta.ampRaw = max(sr(:,bP),[],2);  sta.ampCor = max(sc_(:,bP),[],2);
    sta.ampBld = max(abs(sb(:,bP)),[],2);
    nseR = std(fRaw(:,gd),0,2);  nseC = std(fCor(:,gd),0,2);
    fprintf('[flow] STA over %d spikes (median over ROIs):\n', numel(sIx));
    fprintf('         spike amplitude  %.3g -> %.3g (%+.1f%%)\n', ...
        median(sta.ampRaw), median(sta.ampCor), ...
        100*(median(sta.ampCor)/max(median(sta.ampRaw),eps)-1));
    fprintf('         spike SNR        %.2f -> %.2f (%+.1f%%)\n', ...
        median(sta.ampRaw./nseR), median(sta.ampCor./nseC), ...
        100*(median(sta.ampCor./nseC)/max(median(sta.ampRaw./nseR),eps)-1));
    fprintf('         STA of the REMOVED component: %.3g vs %.3g for the spike itself\n', ...
        median(sta.ampBld), median(sta.ampRaw));
    if median(sta.ampRaw) > 3*median(nseR) && median(sta.ampBld) > 0.15*median(sta.ampRaw)
        warning('flow_compare:staBleed', ...
            ['The removed component has a %.0f%%-of-a-spike deflection in its own STA: ' ...
             'the regressor is partly tracking the AP, not only the vessel.'], ...
            100*median(sta.ampBld)/median(sta.ampRaw));
    end
    % A spike amplitude at or below the noise means the STA measured nothing -- almost
    % always the sign convention (these traces are -F'*mov, i.e. spike-POSITIVE; a movie
    % whose spikes are already positive gives -amplitude here and max() then picks noise).
    % Report that rather than dividing by it and printing a ratio with twenty digits.
    okAmp = median(sta.ampRaw) > 3*median(nseR);
    if ~okAmp
        warning('flow_compare:noSpikeAmp', ...
            ['The spike-triggered average has no measurable positive peak (median %.3g ' ...
             'against noise %.3g). These traces are built as -F''*mov and must be ' ...
             'spike-POSITIVE; check the sign of data.mov_res. STA-based normalisation is ' ...
             'switched off.'], median(sta.ampRaw), median(nseR));
        sta.bapMeas = [];
    else
        sta.bapMeas = sta.ampRaw;
    end
else
    fprintf('[flow] only %d usable spikes; the STA check is skipped\n', numel(sIx));
end

% ---- normalisation for display --------------------------------------------------------
[scale, normLbl] = traceScale(cfg, sta, fRaw, gd);

% ---- 2. where in frequency ------------------------------------------------------------
nW2 = min(c.psdWin, 2^floor(log2(nT/2)));
hop = round(nW2/2);  han = 0.5*(1-cos(2*pi*(0:nW2-1)/(nW2-1)));
st  = 1:hop:(nT-nW2+1);  fvec = (0:floor(nW2/2))*fs/nW2;
Praw = zeros(nR, numel(fvec));  Pcor = Praw;
for b = 1:numel(st)
    ib = st(b) + (0:nW2-1);
    for q = 1:2
        X = dRaw(:,ib);  if q==2, X = dCor(:,ib); end
        X = X - mean(X,2);  X(:,~gd(ib)) = 0;
        Z = abs(fft(X.*han, [], 2)).^2;  Z = Z(:,1:numel(fvec));
        if q==1, Praw = Praw + Z; else, Pcor = Pcor + Z; end
    end
end
Praw = Praw/numel(st);  Pcor = Pcor/numel(st);

cp = struct('Traw',Traw, 'Tcor',Tcor, 'Tbld',Tbld, 'dRaw',dRaw, 'dCor',dCor, 'dBld',dBld, ...
    'remFull',remAll, 'remFast',remFast, 'remSub',remSub, 'onVessel',onVes, ...
    'bvFrac',bvFrac, 'segFrac',segFrac, 'onSeg',onSeg, 'segMask',segMask, ...
    'varBld',var(dBld(:,gd),0,2), 'varRaw',var(dRaw(:,gd),0,2), ...
    'sta',sta, 'scale',scale, 'normLbl',normLbl, ...
    'psdF',fvec, 'psdRaw',Praw, 'psdCor',Pcor, 'good',gd, ...
    'fs',fs, 'sz',[nrow ncol]);
end

% =========================================================================
function [scale, lbl] = traceScale(cfg, sta, fRaw, gd)
% Divide by the bAP measured on the trace's OWN spike-triggered average, so 1.0 = one
% spike and ROIs of different depth are comparable. A floor stops an ROI with no spike
% response from blowing the display up.
c = cfg.compare;
switch c.normMode
    case 'own'
        scale = 1 ./ max(std(fRaw(:,gd),0,2), eps);  lbl = 'own fast-noise SD';
    case 'bap'
        if isempty(sta.bapMeas)
            fprintf(['[flow] no usable bAP, so ''bap'' normalisation is not available; ' ...
                     'falling back to the fast-noise SD\n']);
            scale = 1 ./ max(std(fRaw(:,gd),0,2), eps);  lbl = 'own fast-noise SD';
            return
        end
        b = sta.bapMeas;
        fl = 0.05*median(b,'omitnan');
        nFloor = nnz(b < fl);
        scale = 1 ./ max(b, fl);
        lbl = 'bAP peak (1.0 = one spike)';
        if nFloor > 0
            fprintf('[flow] %d ROI(s) had a bAP below 5%% of the median and were floored\n', nFloor);
        end
    otherwise
        error('cfg.compare.normMode must be ''bap'' or ''own''.');
end
end
