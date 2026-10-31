function bp = flow_bandpass(data, cfg)
% bp = flow_bandpass(data, cfg)
%
% Spatial band-pass of the motion-residual movie: keep only structure on the vessel's
% spatial scale.
%
% WHY SPATIALLY. The three things living in mov_res sit at very different spatial scales:
%     shot noise    independent per pixel        (< 1 um)
%     blood         coherent over ~10 um, gone by 30 um
%     dendritic Vm  coherent over the length constant  (>> 50 um)
% so a purely SPATIAL band-pass isolates the vascular component using no temporal
% assumption at all -- no filter band, no lag search, no clustering. Whatever survives is
% by construction the part of the movie that varies on the vessel's scale.
%
% VOLTAGE FIRST, AND FROM BESIDE THE VESSELS. The dendrite is a thin bright structure, so
% the filter responds to its EDGES and a band-pass alone only partly suppresses Vm.
% Something has to be projected out first -- but WHAT is the question cfg.bandpass.voltRef
% answers, and the two answers are not equivalent:
%
%   'traces' uses data.voltTr, the per-ROI voltage traces. Those are read off the WHOLE
%            footprint, and the footprint OVERLAPS the vessels, so that "voltage" already
%            contains blood. Regressing it out takes vascular signal away before anything
%            has measured it, and the estimate removed is itself part vascular.
%   'ring'   builds the reference where blood cannot be: a ring of footprint just outside
%            the vessel mask, split into connected patches so that two pieces of dendrite
%            at different points on the neurite never get averaged together. See
%            get_ringRef. This is the default.
%
% The ring is drawn around data.bvBW, NOT around the detection-density mask. densBW would
% be the sharper choice -- it is where RBCs were actually seen -- but it is produced by
% flow_segments, three stages downstream, and using it here would make the pipeline
% circular. bvBW is an input, available before anything runs, and covers vessels the
% tracking never found.
%
% THE ACHIEVED BAND IS MEASURED, NOT ASSUMED. A Difference of Gaussians cannot make a
% band narrower than ~3.4:1 whatever sigmas are used, so a request for 6-30 um comes out
% wider than asked. The radial transfer function is computed and its half-power band
% reported, and a warning fires if it is far from the request. Read bp.bandUm, not
% cfg.bandpass.loUm/hiUm, when describing what the filter did.
%
% INPUT   data  .mov_res .pxUm .bvBW .neuronMap
%               .ftprnt   required for voltRef = 'ring'
%               .voltTr   required for voltRef = 'traces'
%         cfg   uses cfg.bandpass and cfg.common.fs
% OUTPUT  bp
%   .S        [nrow x ncol x nT] single, the band-passed movie
%   .mov_nv   [nrow x ncol x nT] single, the input to the filter (voltage-regressed)
%   .volt     what was removed and where: .mode .n .Tr the traces, .W their footprint
%             weights (one column each, in the same order as .Tr), .info the ring map and
%             the patch bookkeeping. .W is what to draw to see what the movie actually
%             lost -- it is carried through the dedup and the cap, so it is the survivors
%             and not the candidates.
%   .bandUm   [lo hi] MEASURED half-power band, um
%   .peakUm   the period the filter passes best
%   .kernel   the DoG kernel, or the Fourier mask, whichever was used
%   .radial   struct with .perUm and .gain, the measured transfer function
%   .varRes .varS   per-pixel variance before and after
%   .stats    struct of the concentration ratios printed below
%
% 2026.08.28, Byung Hun Lee

c = cfg.bandpass;
[nrow, ncol, nT] = size(data.mov_res);
pxUm = data.pxUm;

% ---- remove the neuron's voltage signal ---------------------------------------------
% WHICH VESSEL MASK THE RING IS DRAWN AROUND. In the script this stage came after the
% tracking, so the ring could be drawn around densBW -- the mask of where RBCs were
% actually detected. Here it cannot: densBW is produced by flow_segments, which is three
% stages downstream of this one, and waiting for it would make the pipeline circular.
% data.bvBW is the anatomical vessel mask, is an INPUT rather than a pipeline product, and
% is available before anything runs. It is also the more conservative choice -- it covers
% vessels the tracking never found, and every pixel it covers is one the ring must avoid.
volt = struct('mode', lower(c.voltRef), 'n', 0, 'W', [], 'Tr', [], 'info', []);
switch lower(c.voltRef)

case 'ring'
    assert(~isempty(data.ftprnt), ['cfg.bandpass.voltRef = ''ring'' needs data.ftprnt: ' ...
        'the reference is the part of the footprint that falls in a ring outside the ' ...
        'vessels, and with no footprint there is nothing to intersect. Supply ftprnt, ' ...
        'or set voltRef to ''traces'' or ''none''.']);
    assert(any(data.bvBW(:)), ['cfg.bandpass.voltRef = ''ring'' needs data.bvBW: the ring ' ...
        'is drawn around it. Supply the vessel mask, or set voltRef to ''traces'' or ' ...
        '''none''.']);
    smN = max(1, round(c.ringSmoothMs * cfg.common.fs / 1000));
    [Vr, Wr, rinf] = get_ringRef(data.mov_res, data.ftprnt, data.bvBW, struct( ...
        'ringIn', c.ringIn, 'ringOut', c.ringOut, 'thr', c.ringThr, ...
        'minPx', c.ringMinPx, 'smoothN', smN, 'dedupR', c.ringDedupR, ...
        'maxTrace', c.ringMax, 'chunk', c.chunk));
    fprintf(['[flow] bandpass: ring reference %d-%d px outside bvBW, %d px of ring; ' ...
             '%d patches\n'], c.ringIn, c.ringOut, nnz(rinf.ring), rinf.nPatch);
    fprintf(['         -> %d traces regressed out (%d duplicates, %d over the cap), ' ...
             'smoothed %g ms\n'], size(Vr,1), rinf.nDup, rinf.nCapped, c.ringSmoothMs);
    mov_nv = single(SeeResiduals(data.mov_res, Vr));
    volt.n = size(Vr,1);  volt.W = Wr;  volt.Tr = Vr;  volt.info = rinf;
    % The number to watch: the reference must NOT be taking the vascular signal with it.
    v0 = var(double(data.mov_res), 0, 3);  v1 = var(double(mov_nv), 0, 3);
    if any(data.bvBW(:))
        keptBv = mean(v1(data.bvBW)) / max(mean(v0(data.bvBW)), eps);
        fprintf('         variance kept on bvMask %.1f%%, off it %.1f%%\n', ...
            100*keptBv, 100*mean(v1(~data.bvBW))/max(mean(v0(~data.bvBW)), eps));
        if keptBv < 0.8
            warning('flow_bandpass:ringAteVessel', ...
                ['The ring reference removed %.0f%% of the variance ON THE VESSELS. That ' ...
                 'is the signal this pipeline exists to measure -- the ring is probably ' ...
                 'touching vessel pixels. Raise cfg.bandpass.ringIn.'], 100*(1-keptBv));
        end
    end
    clear v0 v1

case 'traces'
    assert(~isempty(data.voltTr), ['cfg.bandpass.voltRef = ''traces'' but data.voltTr is ' ...
        'empty. Use ''ring'' (which needs only ftprnt and bvBW) or ''none''.']);
    V = data.voltTr;
    % Accept either orientation. tovec(mov)'*tovec(ftprnt) -- the natural way to write
    % it -- comes out [nT x nROI], while everything here works in [nROI x nT], and a
    % silent mismatch would either error deep inside SeeResiduals or, when nROI happens
    % to equal nT, regress out the wrong thing entirely.
    if size(V,2) ~= nT && size(V,1) == nT
        V = V.';
        fprintf('[flow] bandpass: voltTr given as [nT x nTrace]; transposed to [%d x %d]\n', ...
            size(V,1), size(V,2));
    end
    assert(size(V,2) == nT, ['voltTr is %s but the movie has %d frames -- it must be ' ...
        '[nTrace x nFrames] or [nFrames x nTrace].'], mat2str(size(data.voltTr)), nT);
    keepV = all(isfinite(V), 2) & std(V, 0, 2) > 0;
    V = V(keepV, :);
    V = (V - mean(V,2)) ./ std(V,0,2);             % scale-free, better conditioned
    fprintf(['[flow] bandpass: regressing %d of %d voltTr traces out of mov_res. Those ' ...
             'are read off the\n         WHOLE footprint, which overlaps the vessels, so ' ...
             'some vascular signal goes with\n         them -- voltRef = ''ring'' avoids ' ...
             'that.\n'], nnz(keepV), numel(keepV));
    mov_nv = single(SeeResiduals(data.mov_res, V));
    volt.n = size(V,1);  volt.Tr = V;

case 'none'
    fprintf(['[flow] bandpass: no voltage removed; the spatial filter does all the work ' ...
             '(expect Vm at the dendrite edges)\n']);
    mov_nv = single(data.mov_res);

otherwise
    error('cfg.bandpass.voltRef must be ''ring'', ''traces'' or ''none''.');
end

% ---- build the filter ----------------------------------------------------------------
fLoC = 1/c.hiUm;  fHiC = 1/c.loUm;                 % passband in cycles per um
kDoG = [];  Hmask = [];
switch lower(c.method)
    case 'fourier'
        [fx, fy] = meshgrid(((0:ncol-1)-floor(ncol/2))/ncol, ...
                            ((0:nrow-1)-floor(nrow/2))/nrow);
        fr = sqrt(fx.^2 + fy.^2) / pxUm;
        Hband = 1 ./ (1 + (fr/fHiC).^(2*c.order)) .* ...
                1 ./ (1 + (fLoC./max(fr,eps)).^(2*c.order));
        Hmask  = ifftshift(Hband);
        frMeas = fr;  Hmeas = Hband;  kern = Hband;
    case 'dog'
        % Calibrate the two sigmas numerically so the MEASURED band is centred on the
        % request. The analytic sigma = d/2 rule is wrong by ~2.7x for a Gaussian, whose
        % half-power period is ~10.7*sigma rather than 2*sigma.
        rat = 1.2;  a = 0.2;  b = 20;  tgt = sqrt(c.loUm*c.hiUm);  s1 = sqrt(a*b);
        for it = 1:40
            s1 = sqrt(a*b);
            kk = fspecial('gaussian', 2*ceil(4*max(s1*rat,1))+1, s1) - ...
                 fspecial('gaussian', 2*ceil(4*max(s1*rat,1))+1, s1*rat);
            KK = abs(fftshift(fft2(kk, 256, 256)));
            [gx, gy] = meshgrid((-128:127)/256);  gr = sqrt(gx.^2+gy.^2);
            eE = linspace(0,0.5,200);  eC = (eE(1:end-1)+eE(2:end))/2;
            Kr = arrayfun(@(q) mean(KK(gr>=eE(q) & gr<eE(q+1)),'omitnan'), 1:numel(eC));
            [~, ip] = max(Kr);  pk = 1/(eC(ip)/pxUm);
            if pk > tgt, b = s1; else, a = s1; end
        end
        sLo = s1;  sHi = s1*rat;
        kDoG = fspecial('gaussian', 2*ceil(4*sHi)+1, sLo) - ...
               fspecial('gaussian', 2*ceil(4*sHi)+1, sHi);
        [fx, fy] = meshgrid((-128:127)/256);
        frMeas = sqrt(fx.^2+fy.^2) / pxUm;
        Hmeas  = abs(fftshift(fft2(kDoG, 256, 256)));
        kern   = kDoG;
        fprintf('[flow] bandpass: DoG calibrated, sigma %.2f / %.2f px\n', sLo, sHi);
    otherwise
        error('cfg.bandpass.method must be ''dog'' or ''fourier''.');
end

% ---- what band does this filter ACTUALLY pass? ---------------------------------------
eE = linspace(0, max(frMeas(:)), 300);  eC = (eE(1:end-1)+eE(2:end))/2;
Hrad = arrayfun(@(q) mean(Hmeas(frMeas>=eE(q) & frMeas<eE(q+1)), 'omitnan'), 1:numel(eC));
ok   = isfinite(Hrad) & eC > 0;
perUm = 1 ./ eC(ok);  gain = Hrad(ok);
[hPk, ip] = max(gain);  half = gain >= hPk/2;
bandUm = [min(perUm(half)) max(perUm(half))];
fprintf('[flow] bandpass: requested %g-%g um, ACHIEVED %.1f-%.1f um (peak %.1f, ratio %.2f vs %.2f)\n', ...
    c.loUm, c.hiUm, bandUm(1), bandUm(2), perUm(ip), bandUm(2)/bandUm(1), c.hiUm/c.loUm);
if bandUm(1) < 0.7*c.loUm || bandUm(2) > 1.4*c.hiUm
    warning('flow_bandpass:wideBand', ...
        ['The achieved band is much wider than requested -- interpret S as the ' ...
         '%.1f-%.1f um band, not %g-%g um.'], bandUm(1), bandUm(2), c.loUm, c.hiUm);
end

% ---- apply it -------------------------------------------------------------------------
S = zeros(nrow, ncol, nT, 'single');
for c0 = 1:c.chunk:nT
    ci = c0:min(c0+c.chunk-1, nT);
    if strcmpi(c.method, 'fourier')
        S(:,:,ci) = single(real(ifft2(fft2(mov_nv(:,:,ci)) .* Hmask)));
    else
        S(:,:,ci) = imfilter(mov_nv(:,:,ci), kDoG, 'replicate');
    end
end

% ---- does S concentrate on the vessels? ----------------------------------------------
varRes = var(double(mov_nv), 0, 3);
varS   = var(double(S), 0, 3);
bvp    = data.bvBW(:) > 0;   neu = data.neuronMap(:) > 0;
% a ratio needs pixels on BOTH sides; an absent mask reports NaN rather than a number
% that looks meaningful
rat = @(v, m) merge(any(m) && ~all(m), mean(v(m))/max(mean(v(~m)), eps), NaN);
st = struct();
st.bvRatioIn  = rat(varRes, bvp);   st.bvRatioS  = rat(varS, bvp);
st.neuRatioIn = rat(varRes, neu);   st.neuRatioS = rat(varS, neu);
am = true(size(varS));
if isfield(data,'anaMask') && isequal(size(data.anaMask), size(varS)), am = data.anaMask; end
st.keptPct    = 100*mean(varS(am))/max(mean(varRes(am)), eps);   % inside the mask
st.anaMask    = am;
if isfinite(st.bvRatioIn)
    fprintf('       variance in bvMask / outside:  input %.2f -> S %.2f\n', st.bvRatioIn, st.bvRatioS);
else
    fprintf('       (no bvMask supplied, so the vessel concentration ratio is not computed)\n');
end
if isfinite(st.neuRatioIn)
    fprintf('       variance on neuron / off:      input %.2f -> S %.2f\n', st.neuRatioIn, st.neuRatioS);
end
fprintf('       S keeps %.1f%% of the input variance inside the footprint mask (%d px)\n', ...
    st.keptPct, nnz(am));

bp = struct('S', S, 'mov_nv', mov_nv, 'bandUm', bandUm, 'peakUm', perUm(ip), ...
            'kernel', kern, 'radial', struct('perUm', perUm, 'gain', gain), ...
            'varRes', varRes, 'varS', varS, 'stats', st, 'method', c.method, ...
            'volt', volt);
end

function v = merge(cond, a, b)
if cond, v = a; else, v = b; end
end
