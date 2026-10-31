function [Tr, W, info] = get_ringRef(mov, ftprnt, vesselBW, opt)
% [Tr, W, info] = get_ringRef(mov, ftprnt, vesselBW, opt)
%
% Voltage reference traces taken from BESIDE the vessels, where blood cannot be.
%
% THE PROBLEM THIS SOLVES. The obvious voltage reference is the per-ROI trace read off the
% whole footprint. But the footprint overlaps the vessels, so that "voltage" already has
% blood in it: regressing it out of the movie takes some of the vascular signal away
% before anything has measured it, and leaves a voltage estimate that is itself part
% vascular.
%
% Here the reference is built in a RING just outside each vessel, intersected with the
% neuron's footprint. Those pixels are dendrite adjacent to the vessel -- close enough to
% carry the same membrane potential, far enough to carry none of that vessel's RBC train.
%
% ONE ROI'S OVERLAP WITH ONE RING IS USUALLY SEVERAL SEPARATE PATCHES, and they must not
% be averaged together. A vessel can be far longer than the stretch of dendrite carrying
% a given ROI's signal, so the ring around it clips that footprint in two, three, four
% disconnected places -- on opposite sides of the vessel, or at either end. Those patches
% sit at different points along the neurite and so are not at the same potential; one mean
% over all of them is a reference that matches none of them. Each connected patch above
% opt.minPx therefore becomes its OWN trace, and patches below it are dropped rather than
% merged (a handful of pixels is shot noise, and regressing shot noise out of a movie
% removes variance everywhere for nothing).
%
% THE SMOOTHING IS PART OF THE ARGUMENT. A ring trace is dendritic Vm plus shot noise at
% every frequency. opt.smoothN keeps the part that is actually voltage. Note what it does
% NOT do: a 5 ms boxcar is -3 dB near 90 Hz, so it does not confine the reference to the
% sub-40 Hz voltage band -- it still passes power through the lower half of the RBC range.
% That is acceptable ONLY because the ring excludes every vessel pixel, so there is no RBC
% train in the trace for it to pass. Widen the ring, or drop the exclusion, and this has
% to grow a lot.
%
% INPUT
%   mov       [nrow x ncol x nT]  the movie the traces are read from (e.g. mov_res)
%   ftprnt    [nrow x ncol x nROI] per-ROI footprints; [] is an error, there is nothing
%             to intersect the ring with
%   vesselBW  logical [nrow x ncol]  what the ring is drawn AROUND, and what is excluded
%             from it. Pass the UNION of every vessel you know about -- a ring around
%             vessel A that lands on vessel B would put blood straight back in.
%   opt   struct, all optional:
%     .ringIn    px, gap between the mask and the inner edge of the ring   (default 1)
%                A mask that is a threshold on a smoothed map has a soft boundary: the
%                pixel just outside it still carries vessel signal.
%     .ringOut   px, outer edge of the ring                                (default 6)
%     .thr       fraction of an ROI's own peak below which its footprint
%                counts as absent                                         (default 0.2)
%     .minPx     px, smallest connected patch that becomes a trace         (default 50)
%     .smoothN   frames of moving average applied to each trace            (default 5)
%     .dedupR    |corr| above which two traces count as the same one       (default 0.95)
%                Rings around neighbouring vessels overlap the same piece of dendrite and
%                would otherwise contribute it several times over.
%     .maxTrace  hard cap on how many traces are returned                  (default 40)
%                Every trace is a degree of freedom taken out of the movie.
%     .exclude   extra logical mask to keep out of every ring              (default none)
%     .chunk     frames per pass when reading the movie                    (default 500)
%
% OUTPUT
%   Tr    [nTrace x nT]  z-scored reference traces, ready for SeeResiduals
%   W     [nPix x nTrace] the footprint weights each trace came from, normalised to sum 1
%         and in the SAME order as Tr, so a picture of W is a picture of what was removed
%   info  .ring  logical map of the union of the rings actually used
%         .nPatch  how many patches were found before dedup and the cap
%         .nDup .nCapped  how many were dropped by each
%         .roiOf .patchOf  which ROI and which patch of it each trace came from
%
% 2026.08.30, Byung Hun Lee

if nargin < 4, opt = struct; end
ringIn   = pick_(opt, 'ringIn',   1);
ringOut  = pick_(opt, 'ringOut',  6);
thr      = pick_(opt, 'thr',      0.2);
minPx    = pick_(opt, 'minPx',    50);
smoothN  = pick_(opt, 'smoothN',  5);
dedupR   = pick_(opt, 'dedupR',   0.95);
maxTrace = pick_(opt, 'maxTrace', 40);
chunk    = pick_(opt, 'chunk',    500);
excl     = pick_(opt, 'exclude',  []);

[nrow, ncol, nT] = size(mov);
npx = nrow * ncol;
assert(~isempty(ftprnt), ['get_ringRef needs per-ROI footprints: there is nothing for ' ...
    'the ring to intersect with. Pass data.ftprnt, or use a different voltage reference.']);
assert(isequal(size(ftprnt,1), nrow) && isequal(size(ftprnt,2), ncol), ...
    'ftprnt is %s but the movie frame is %s.', mat2str(size(ftprnt)), mat2str([nrow ncol]));
assert(isequal(size(vesselBW), [nrow ncol]), 'vesselBW is %s but the frame is %s.', ...
    mat2str(size(vesselBW)), mat2str([nrow ncol]));
assert(any(vesselBW(:)), ['vesselBW selects no pixel, so there is nothing to draw a ring ' ...
    'around.']);
keepOut = logical(vesselBW);
if ~isempty(excl), keepOut = keepOut | logical(excl); end

% ---- the ring, and the patches of footprint inside it -------------------------------
if ringIn > 0, inner = imdilate(vesselBW, strel('disk', ringIn)); else, inner = vesselBW; end
ring = imdilate(vesselBW, strel('disk', ringOut)) & ~inner & ~keepOut;
assert(any(ring(:)), ['The ring is empty: ringOut = %d px does not clear ringIn = %d px ' ...
    'anywhere outside the excluded masks.'], ringOut, ringIn);

W = zeros(npx, 0);  roiOf = [];  patchOf = [];
for j = 1:size(ftprnt, 3)
    fj = double(ftprnt(:,:,j));
    mx = max(fj(:));
    if ~isfinite(mx) || mx <= 0, continue, end
    ov = (fj >= thr * mx) & ring;
    if ~any(ov(:)), continue, end
    CCo = bwconncomp(ov);
    for p = 1:CCo.NumObjects
        pxo = CCo.PixelIdxList{p};
        if numel(pxo) < minPx, continue, end
        w = zeros(nrow, ncol);
        w(pxo) = fj(pxo);                       % keep the footprint's own weighting
        W(:, end+1)  = w(:) / sum(w(:));        %#ok<AGROW>
        roiOf(end+1)   = j;                     %#ok<AGROW>
        patchOf(end+1) = p;                     %#ok<AGROW>
    end
end
nPatch = size(W, 2);
assert(nPatch > 0, ['No ring-and-footprint overlap held a connected patch of at least ' ...
    'minPx = %d px. Widen ringOut, lower thr or minPx, or check that the footprints ' ...
    'actually reach the vessels in this field.'], minPx);

% ---- the traces: one chunked pass, all patches at once -------------------------------
% Chunked because a full [nPix x nT] copy of the movie is several GB. The weights are
% sparse, so the product is cheap once a chunk is in memory.
Tr = zeros(nPatch, nT);
for c0 = 1:chunk:nT
    ci = c0:min(c0+chunk-1, nT);
    Tr(:, ci) = W' * double(reshape(mov(:,:,ci), npx, numel(ci)));
end

if smoothN > 1, Tr = movmean(Tr, smoothN, 2); end
ok = all(isfinite(Tr), 2) & std(Tr, 0, 2) > 0;
Tr = Tr(ok,:);  W = W(:,ok);  roiOf = roiOf(ok);  patchOf = patchOf(ok);
Tr = (Tr - mean(Tr,2)) ./ std(Tr,0,2);

% ---- drop duplicates, then cap. W is carried through every step, so what comes back
% is exactly the set that will be regressed out -- a picture of the candidates rather
% than of the survivors would be a picture of something the movie never saw.
keep = true(1, size(Tr,1));
for a = 2:size(Tr,1)
    prev = find(keep(1:a-1));
    if isempty(prev), continue, end
    if any(abs(corr(Tr(a,:)', Tr(prev,:)')) > dedupR), keep(a) = false; end
end
nDup = nnz(~keep);
Tr = Tr(keep,:);  W = W(:,keep);  roiOf = roiOf(keep);  patchOf = patchOf(keep);

nCapped = 0;
if size(Tr,1) > maxTrace
    [~, ord] = sort(var(Tr,0,2), 'descend');    % the ones actually carrying signal
    ord = sort(ord(1:maxTrace));
    nCapped = size(Tr,1) - maxTrace;
    Tr = Tr(ord,:);  W = W(:,ord);  roiOf = roiOf(ord);  patchOf = patchOf(ord);
end

info = struct('ring', ring, 'nPatch', nPatch, 'nDup', nDup, 'nCapped', nCapped, ...
    'roiOf', roiOf, 'patchOf', patchOf, 'ringIn', ringIn, 'ringOut', ringOut, ...
    'minPx', minPx, 'smoothN', smoothN);
end

% -------------------------------------------------------------------------------------
function v = pick_(s, f, dflt)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
