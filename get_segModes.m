function [U, tr, sc, nm] = get_segModes(Y, fs, band, method, opt)
% [U, tr, sc, nm] = get_segModes(Y, fs, band, method, opt)
%
% Decompose one vessel segment into candidate oscillatory modes, three ways.
%
% Y is the segment's pixels over time. What comes back is a set of MODES: a spatial
% weight vector, a time course, and an oscillation score for each, already filtered to
% the ones that look like an oscillation in `band` and ordered best first. The caller
% turns each time course into a cos/sin pair.
%
% WHY THREE METHODS. They differ in what they assume about the segment.
%
%   'cpca'  COMPLEX PCA. Band-limit, take the analytic signal, SVD. A travelling wave is
%           rank 1 in that domain instead of rank 2, so one mode holds the whole wave --
%           |U| the envelope along the vessel, angle(U) the propagation. Splits the
%           segment by VARIANCE, which is the right criterion when the segment holds one
%           dominant source and the rest is noise.
%
%   'ica'   PCA-reduce, then FastICA on the component TIME COURSES. Variance is the wrong
%           criterion when a connected component holds TWO capillaries of similar
%           amplitude: the leading PCs are then mixtures of both, and neither is a vessel.
%           Independence separates them, and two RBC trains at different speeds are close
%           to independent. What it costs: ICA has no ordering, no orthogonality, and no
%           complex form here -- so a travelling wave can come back as two real ICs 90 deg
%           apart, i.e. one physical source spread over two modes. The caller should check
%           the design's collinearity.
%
%   'disk'  NO decomposition. Tile the segment with overlapping disks and take each disk's
%           mean. Every mode is then a local measurement with an obvious spatial meaning
%           and nothing is assumed about how many sources the segment holds. The price is
%           that neighbouring disks overlap by construction, so the modes are correlated
%           and the individual coefficients are not readable -- fair when the question is
%           how much the design REMOVES, not what each row means.
%
% TEMPORAL, NOT SPATIAL ICA. What we want separated is two RBC trains running through one
% component, and what makes them distinct is their time course. Spatial ICA would ask the
% PIXELS to be independent, which for two vessels a few pixels apart they are not.
%
% THE BAND IS APPLIED BEFORE ICA, ALWAYS. ICA maximises non-gaussianity, and over the
% full band the most non-gaussian thing in a segment is the shot-noise tail, not the RBC
% train: unfiltered, the leading IC comes back as noise every time.
%
% INPUT
%   Y       [nPix x nT]  the segment's traces, mean-removed by the caller
%   fs      Hz
%   band    [lo hi] Hz, where the oscillation is expected
%   method  'cpca' | 'ica' | 'disk'
%   opt     struct:
%     .scoreMin   keep a mode only if get_oscScore reaches this      (required)
%     .goodMask   [1 x nT] logical, frames to fit on                 (required)
%     .segId      segment number, for the names                      (required)
%     .maxMode    at most this many modes returned            (default 8)
%     .binPix     px per mode scanned: nScan = ceil(nPix/binPix)     (default 12)
%     .minPxPCA   below this many px, 'cpca' returns the segment mean (default 20)
%     .pcFloor    'ica' needs at least this many PCs to run          (default 3)
%     .pxIdx      [nPix x 1] linear indices, REQUIRED for 'disk'
%     .sz         [nrow ncol], REQUIRED for 'disk'
%     .spacing    px, MINIMUM separation between disk centres        (default 6)
%                 Centres are placed by farthest-point sampling and stop when the next
%                 one would fall closer than this, so it is what actually sets the disk
%                 count -- opt.maxMode is only a ceiling on top of it.
%     .radius     px, disk radius                                    (default 3)
%     .minPxDisk  px of segment inside a disk to use it              (default 7)
%
% OUTPUT
%   U    [nPix x nMode]  spatial weights (complex for 'cpca')
%   tr   [nMode x nT]    time courses, real, defined on ALL frames
%   sc   [nMode x 1]     oscillation score of each
%   nm   {1 x nMode}     a name for each
%   All four are empty when nothing reached opt.scoreMin.
%
% 2026.08.30, Byung Hun Lee

U = [];  tr = [];  sc = [];  nm = {};
if nargin < 5, opt = struct; end
gd  = opt.goodMask(:)';
sgi = opt.segId;
scMin   = opt.scoreMin;
maxMode = pick_(opt, 'maxMode', 8);
binPix  = pick_(opt, 'binPix', 12);
nPx     = size(Y, 1);

switch lower(method)

case 'cpca'
    minPxPCA = pick_(opt, 'minPxPCA', 20);
    if nPx < minPxPCA
        % No room for a phase gradient: the leading mode would be the mean with extra
        % noise, so use the mean and say so.
        tr = mean(Y, 1);  U = ones(nPx,1)/nPx;
        sc = get_oscScore(tr, fs, band);
        nm = {sprintf('seg%d own trace (small)', sgi)};
        if sc < scMin, U = []; tr = []; sc = []; nm = {}; end
        return
    end
    nScan = clampScan_(nPx, binPix, size(Y,1), nnz(gd));
    if nScan < 2, return, end
    [Us, pcTr, varPC] = get_cpca(Y, fs, band, nScan, gd);
    scAll = get_oscScore(real(pcTr), fs, band);
    keep  = find(scAll >= scMin);
    if isempty(keep), return, end
    [~, ord] = sort(scAll(keep) .* varPC(keep), 'descend');   % oscillatory AND strong
    keep = keep(ord);
    if numel(keep) > maxMode, keep = keep(1:maxMode); end
    U  = Us(:, keep);
    tr = real(pcTr(keep, :));
    sc = scAll(keep);
    nm = arrayfun(@(q) sprintf('seg%d cpca%d', sgi, q), keep(:)', 'uni', 0);

case 'ica'
    pcFloor = pick_(opt, 'pcFloor', 3);
    nScan = clampScan_(nPx, binPix, size(Y,1), nnz(gd));
    if nScan < pcFloor, return, end

    Yb = get_bandpass(Y, fs, band);
    Yb(~isfinite(Yb)) = 0;
    Yb = Yb - mean(Yb(:, gd), 2);
    [Us, ~, ~] = svds(Yb(:, gd), nScan);
    pcTr = Us' * Yb;                                    % nScan x nT, ALL frames

    nIC  = min(maxMode, nScan);
    Wun  = [];  Amix = [];
    if ~isempty(which('fastica'))
        try
            % THREE outputs, not two. FastICA overloads its return list by nargout:
            % [icasig, A, W] with three, but [A, W] with two -- asking for two and
            % naming the first one icasig silently hands back the mixing matrix.
            [~, Amix, Wun] = fastica(pcTr(:, gd), 'numOfIC', nIC, ...
                'verbose', 'off', 'displayMode', 'off');
        catch ME
            fprintf('   seg %3d: fastica failed (%s) -> rica\n', sgi, ME.message);
            Wun = [];
        end
    end
    if isempty(Wun)
        if isempty(which('rica')), return, end
        mdl  = rica(pcTr(:, gd)', nIC, 'IterationLimit', 500, 'VerbosityLevel', 0);
        Wun  = mdl.TransformWeights';                   % nIC x nScan
        Amix = pinv(Wun);
    end
    if isempty(Wun) || size(Wun,1) < 1, return, end

    % Wun was estimated on the good frames only; apply it to ALL of them. A regressor
    % with a gap in it cannot be regressed out of the movie, so the components have to be
    % defined everywhere -- and the unmixing is a fixed linear map, so extending it is
    % just a matrix product, not a re-fit.
    icaAll = Wun * pcTr;                                % nIC x nT
    Uall   = Us * Amix;                                 % nPix x nIC

    scAll = get_oscScore(icaAll, fs, band);
    keep  = find(scAll >= scMin);
    if isempty(keep), return, end
    % ICA returns its components in an ARBITRARY order -- there is no eigenvalue to sort
    % by -- so ranking on score alone would put a clean but negligible component first.
    % Weight by how much of the segment each IC actually accounts for, the same way the
    % PCA path weights by explained variance.
    pwr = sum(abs(Uall).^2, 1)' .* var(icaAll, 0, 2);
    pwr = pwr / max(sum(pwr), eps);
    [~, ord] = sort(scAll(keep) .* pwr(keep), 'descend');
    keep = keep(ord);
    if numel(keep) > maxMode, keep = keep(1:maxMode); end
    U  = Uall(:, keep);
    tr = icaAll(keep, :);
    sc = scAll(keep);
    nm = arrayfun(@(q) sprintf('seg%d ic%d', sgi, q), keep(:)', 'uni', 0);

case 'disk'
    assert(isfield(opt,'pxIdx') && isfield(opt,'sz'), ...
        'method ''disk'' needs opt.pxIdx and opt.sz.');
    spacing  = pick_(opt, 'spacing', 6);
    radius   = pick_(opt, 'radius', 3);
    minPxDsk = pick_(opt, 'minPxDisk', 7);
    [rr, cc] = ind2sub(opt.sz, opt.pxIdx(:));

    % CENTRES ARE CHOSEN BY FARTHEST-POINT SAMPLING ON THE SEGMENT'S OWN PIXELS, not by
    % laying a grid over its bounding box. Three things went wrong with the grid:
    %
    %   - a vessel is rarely axis-aligned, so most grid points fall OUTSIDE the segment
    %     and the survivors are wherever the shape happens to line up with the lattice;
    %   - ranking the survivors by how many pixels they contain puts every disk in the
    %     fattest part, because a disk on a tapering end is only half full;
    %   - taking the top-N by SCORE after that ignores position entirely, so a segment
    %     could end up with all of its disks piled on one blob and none along its length.
    %
    % Farthest-point sampling fixes all three at once: seed at the pixel nearest the
    % centroid, then repeatedly take the segment pixel that is farthest from everything
    % already chosen. Every new centre lands in the largest remaining gap, so the centres
    % spread over the segment BY CONSTRUCTION and follow whatever shape it actually has.
    % It stops when the next centre would sit closer than `spacing` to an existing one,
    % which is what makes the disk count scale with the segment instead of being fixed.
    d2near = inf(nPx, 1);
    ctr    = [mean(rr), mean(cc)];
    [~, s0] = min((rr - ctr(1)).^2 + (cc - ctr(2)).^2);   % deterministic seed
    sel = zeros(1, 0);  capBound = true;
    for q = 1:maxMode
        if q == 1
            pick = s0;
        else
            [dmax, pick] = max(d2near);
            if dmax < spacing^2, capBound = false; break, end   % spacing stopped it: good
        end
        sel(end+1) = pick;                        %#ok<AGROW>
        d2near = min(d2near, (rr - rr(pick)).^2 + (cc - cc(pick)).^2);
    end
    if capBound
        % maxMode ran out before `spacing` did, so part of this segment has no disk on it
        % at all. Said out loud, because a bunched tiling looks like a property of the
        % segment when it is really the budget.
        fprintf(['   seg %3d: disk count hit the ceiling (%d) before the %g px spacing ' ...
                 'did -- part of the segment is unsampled; raise diskMaxPer\n'], ...
            sgi, maxMode, spacing);
    end

    loc = zeros(nPx, 0);  keptSel = zeros(1,0);
    for q = 1:numel(sel)
        in = (rr - rr(sel(q))).^2 + (cc - cc(sel(q))).^2 <= radius^2;
        if nnz(in) < minPxDsk, continue, end
        loc(:, end+1)   = in / nnz(in);           %#ok<AGROW>
        keptSel(end+1)  = sel(q);                 %#ok<AGROW>
    end
    if isempty(loc), return, end

    trAll = loc' * Y;
    scAll = get_oscScore(trAll, fs, band);
    keep  = find(scAll >= scMin);
    if isempty(keep), return, end
    % Kept IN PLACE ORDER, not re-sorted by score. The disks are a tiling, and numbering
    % them by rank would scramble the one thing that makes the tiling readable. Whatever
    % the score filter removes leaves a visible hole, which is the honest picture: a gap
    % means no oscillation was found there, not that the disk was never tried.
    if numel(keep) < numel(scAll)
        fprintf('   seg %3d: %d of %d disks below scoreMin, leaving gaps in the tiling\n', ...
            sgi, numel(scAll) - numel(keep), numel(scAll));
    end
    U  = loc(:, keep);
    tr = trAll(keep, :);
    sc = scAll(keep);
    nm = arrayfun(@(q) sprintf('seg%d disk%d', sgi, q), 1:numel(keep), 'uni', 0);

otherwise
    error('method must be ''cpca'', ''ica'' or ''disk''.');
end
end

% -------------------------------------------------------------------------------------
function v = pick_(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end

function n = clampScan_(nPx, binPix, nRow, nGood)
% How many modes to LOOK AT scales with the segment: a longer one spans more independent
% resolution elements and can hold more distinct modes. Scanning a fixed number would
% under-search a big segment and waste time on a small one.
n = max(2, ceil(nPx / binPix));
n = min(n, min(nRow, nGood) - 1);
end
