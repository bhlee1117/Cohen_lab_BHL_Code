function [out, C, info] = SeeResidualsTrace(Y, refSig, varargin)
% [out, C] = SeeResidualsTrace(Y, refSig)
% [out, C] = SeeResidualsTrace(Y, refSig, remOffset)
% [out, C, info] = SeeResidualsTrace(Y, refSig, remOffset, 'fitFrames', m, 'omitnan', true)
%
% SeeResiduals for TRACES: project reference signals out of [nTrace x nT] data instead of
% out of a movie. Same model, same meaning of remOffset, same coefficients -- it just
% never reshapes anything, so it works on traces that were never a movie (ROI traces,
% voltage traces, a single dendrite trace) without packing them into a fake [1 x N x T].
%
% WHAT IT DOES
%   Y(i,:)  ~  C(i,1)*1 + C(i,2)*ref1 + C(i,3)*ref2 + ...
%   out     =  Y - C*I,  the part of each trace that the references cannot explain
%
% THREE THINGS DIFFER FROM SeeResiduals, all because traces are not movies:
%
%  1. NaNs. A movie has none; a trace usually does -- a bad frame, a gap, a blue-light
%     window. SeeResiduals would turn one NaN into a NaN coefficient and hence an
%     all-NaN output trace. Here the fit is taken over the frames that are finite, the
%     coefficients are applied everywhere, and NaNs are put back where the input had
%     them. Set 'omitnan' false to get the strict behaviour.
%
%  2. Frames fitted on need not be frames corrected. 'fitFrames' takes a logical mask:
%     the coefficients come from those frames only and are then applied to all of them.
%     This is the pattern the flow pipeline uses everywhere -- fit on the good tracked
%     window so spikes cannot bias the coefficients, subtract across the whole record --
%     and doing it by hand means slicing, fitting, and re-expanding every time.
%
%  3. The solve. SeeResiduals forms C = Y*I'*inv(I*I') -- an explicit inverse of the
%     normal equations, which squares the condition number and warns or returns garbage
%     when two references are nearly collinear (two ROIs overlapping, a cos and its own
%     derivative). Here it is a QR least-squares solve, mrdivide. For a well-conditioned
%     design the two agree to round-off; for a badly conditioned one this is the one that
%     is still right, and info.rank says when that mattered.
%
% ORIENTATION. Rows are signals and columns are time, as in SeeResiduals. refSig is
% transposed automatically when it arrives as [nT x N]; so is Y, when its orientation is
% unambiguous. It is NOT guessed when the array is square -- that would silently
% transpose a genuine nT x nT case -- and an error says so.
%
% INPUT
%   Y         [nTrace x nT]  the traces to correct (or [nT x nTrace])
%   refSig    [nRef x nT]    the references to remove (or [nT x nRef]).
%             refSig == 0 means "subtract the mean only", exactly as in SeeResiduals.
%   remOffset scalar, default 1. Prepend a row of ones so a constant offset per trace is
%             fitted separately, and centre refSig first so that row means what it says.
%   'fitFrames'  logical [1 x nT], frames the coefficients are estimated from (default all)
%   'omitnan'    fit around NaNs and preserve them in the output (default true)
%
% OUTPUT
%   out   same size AND orientation as Y
%   C     [nTrace x (1+nRef)] coefficients, column 1 the offset when remOffset is on.
%         The analogue of SeeResiduals' scaleImgs, minus the reshape into images.
%   info  .used     [1 x nT] frames the fit actually used (fitFrames minus non-finite)
%         .rank     rank of the design; below size(I,1) means references were collinear
%         .cond     its condition number. Watch this one, not rank: two references that
%                   are collinear to within a whisker still come back FULL rank, while
%                   the coefficients blow up to 1e11 and cancel. Above ~1e8 the residual
%                   is still right and C is meaningless. It warns.
%         .nDropped how many frames were excluded per trace, max over traces
%         .I        the design actually used, [(1+nRef) x nT]
%
% EXAMPLE
%   % take a vessel oscillator out of every ROI trace, fitting only on quiet frames
%   [clean, C] = SeeResidualsTrace(roiTraces, [oscCos; oscSin], 1, 'fitFrames', ~badfr);
%
% 2026.08.30, Byung Hun Lee

% ------------------------------------------------------------------ arguments ---------
remOffset = 1;
if ~isempty(varargin) && (isnumeric(varargin{1}) || islogical(varargin{1})) ...
        && isscalar(varargin{1})
    remOffset = varargin{1};  varargin(1) = [];
end
p = struct('fitframes', [], 'omitnan', true);
for k = 1:2:numel(varargin)
    f = lower(varargin{k});
    assert(isfield(p,f), 'Unknown option ''%s''.', varargin{k});
    p.(f) = varargin{k+1};
end
assert(isscalar(remOffset), 'remOffset must be a scalar (zero or nonzero).');
assert(~isempty(Y) && ismatrix(Y), 'Y must be a non-empty [nTrace x nT] matrix.');

Yin = Y;
flipped = false;

% ---- the "subtract the mean only" case, kept identical to SeeResiduals ---------------
if isscalar(refSig) && refSig == 0
    if size(Y,2) == 1 && size(Y,1) > 1, Y = Y.';  flipped = true; end
    C   = mean(Y, 2, 'omitnan');
    out = Y - C;
    if p.omitnan, out(~isfinite(Y)) = NaN; end
    info = struct('used', all(isfinite(Y),1), 'rank', 1, 'nDropped', 0, ...
                  'I', ones(1, size(Y,2)));
    if flipped, out = out.'; end
    return
end

% ---- orient refSig, then Y, against it -----------------------------------------------
assert(ismatrix(refSig) && ~isempty(refSig), 'refSig must be a matrix, or the scalar 0.');
[rA, rB] = size(refSig);
if size(Y,2) == rB && size(Y,2) ~= rA
    nT = rB;                                   % refSig is already [nRef x nT]
elseif size(Y,2) == rA && size(Y,2) ~= rB
    refSig = refSig.';  nT = rA;
elseif size(Y,1) == rB && size(Y,1) ~= rA
    Y = Y.';  flipped = true;  nT = rB;
elseif size(Y,1) == rA && size(Y,1) ~= rB
    refSig = refSig.';  Y = Y.';  flipped = true;  nT = rA;
elseif rA == rB
    % square refSig: the only consistent reading is that time is whichever dimension Y
    % already agrees with, and if Y agrees with both there is nothing to go on.
    assert(size(Y,2) == rA || size(Y,1) == rA, ...
        ['refSig is %s and Y is %s: they share no dimension, so neither can be time.'], ...
        mat2str(size(refSig)), mat2str(size(Y)));
    assert(~(size(Y,1) == rA && size(Y,2) == rA), ...
        ['refSig is square (%dx%d) and so is Y, so which dimension is time is ' ...
         'ambiguous. Pass them as [nRef x nT] and [nTrace x nT] explicitly.'], rA, rB);
    if size(Y,1) == rA, Y = Y.';  flipped = true; end
    nT = rA;
else
    error(['refSig is %s and Y is %s: they share no dimension that could be time. ' ...
           'Expected [nRef x nT] and [nTrace x nT].'], ...
           mat2str(size(refSig)), mat2str(size(Y)));
end
[nY, ~] = size(Y);
nRef = size(refSig, 1);

% ---- the design ----------------------------------------------------------------------
if remOffset
    dRef = refSig - mean(refSig, 2, 'omitnan');   % centre, so the ones row IS the offset
    I = [ones(1, nT); dRef];
else
    I = refSig;
end

fitM = true(1, nT);
if ~isempty(p.fitframes)
    fitM = logical(p.fitframes(:)');
    assert(numel(fitM) == nT, 'fitFrames has %d elements but there are %d frames.', ...
        numel(fitM), nT);
end
okI = all(isfinite(I), 1) & fitM;      % a non-finite reference frame is unusable for all
assert(any(okI), ['No frame is usable: fitFrames and the finite part of refSig do not ' ...
    'overlap anywhere.']);

% ---- solve -----------------------------------------------------------------------------
% Grouped by NaN pattern rather than looped per trace: traces from the same source
% usually share their bad frames exactly, so the common case is one solve for everything,
% and the loop only costs something when the patterns genuinely differ.
C = zeros(nY, size(I,1));
finY = isfinite(Y);
if ~p.omitnan, finY = true(size(Y)); end
sameAll = all(finY(:)) || all(all(finY == finY(1,:), 1));
nDropped = 0;
if sameAll
    g = okI & (finY(1,:) | ~p.omitnan);
    assert(nnz(g) > size(I,1), ...
        ['Only %d usable frames for %d design rows -- the fit is underdetermined. ' ...
         'Widen fitFrames or drop references.'], nnz(g), size(I,1));
    Yg = Y(:, g);
    if ~p.omitnan, Yg(~isfinite(Yg)) = 0; end
    C = Yg / I(:, g);
    nDropped = nT - nnz(g);
    rk = rank(I(:, g));  cn = cond(I(:, g));
else
    rk = rank(I(:, okI));  cn = cond(I(:, okI));
    for i = 1:nY
        g = okI & finY(i,:);
        if nnz(g) <= size(I,1)
            C(i,:) = NaN;  nDropped = max(nDropped, nT - nnz(g));
            continue
        end
        C(i,:) = Y(i, g) / I(:, g);
        nDropped = max(nDropped, nT - nnz(g));
    end
end
% Rank alone is not enough. Two references that are collinear to within a whisker of
% round-off still come back full rank at rank()'s default tolerance, while the solve
% happily returns coefficients of order 1e11 that cancel against each other. The residual
% is fine either way -- that is what least squares guarantees -- but C is meaningless, and
% the condition number is the thing that says so.
if rk < size(I,1)
    warning('SeeResidualsTrace:rankDeficient', ...
        ['The design has rank %d for %d rows: at least two references are collinear. ' ...
         'The residual is still well defined, but the coefficients in C are not.'], ...
        rk, size(I,1));
elseif cn > 1e8
    warning('SeeResidualsTrace:illConditioned', ...
        ['The design condition number is %.3g: two references are very nearly ' ...
         'collinear. The residual is still well defined, but the coefficients in C are ' ...
         'large and opposed and mean nothing individually (max |C| = %.3g).'], ...
        cn, max(abs(C(:))));
end

out = Y - C*I;
if p.omitnan, out(~finY) = NaN; end

info = struct('used', okI, 'rank', rk, 'cond', cn, 'nDropped', nDropped, 'I', I);
if flipped, out = out.'; end
if ~isequal(size(out), size(Yin))         % should not happen; a loud check beats a silent
    error('SeeResidualsTrace: output is %s but the input was %s.', ...
        mat2str(size(out)), mat2str(size(Yin)));
end
end
