function data = flow_checkData(data)
% data = flow_checkData(data)
%
% Validate and normalise the recording struct ONCE, at the entrance, so a mismatch is
% reported here with the field name rather than surfacing as "Arrays have incompatible
% sizes" four stages later inside a mask expression.
%
% Fills in what is missing, collapses a 3-D mask to 2-D, fixes trace orientation, and
% refuses anything whose size cannot be reconciled with the movie.
%
% 2026.08.28, Byung Hun Lee

assert(isstruct(data) && isfield(data,'mov_res'), 'data.mov_res is required.');
assert(ndims(data.mov_res) == 3, 'data.mov_res must be [nrow x ncol x nFrames], got %s.', ...
    mat2str(size(data.mov_res)));
[nrow, ncol, nT] = size(data.mov_res);
sz = [nrow ncol];

% ---- scalars -------------------------------------------------------------------------
data = fill(data, 'fs',   1000);
data = fill(data, 'pxUm', 1);
data = fill(data, 'tag',  '');
assert(isscalar(data.fs)   && data.fs   > 0, 'data.fs must be a positive scalar.');
assert(isscalar(data.pxUm) && data.pxUm > 0, 'data.pxUm must be a positive scalar.');

% ---- 2-D masks and maps ---------------------------------------------------------------
% A 3-D stack here is almost always the footprint set passed where its projection was
% wanted, so collapse it rather than failing -- but say so, because it is worth knowing.
for f = {'bvBW','neuronMap'}
    nm = f{1};
    dflt = zeros(sz);  if strcmp(nm,'bvBW'), dflt = false(sz); end
    data = fill(data, nm, dflt);
    v = data.(nm);
    if isempty(v), data.(nm) = dflt; continue, end
    if ndims(v) == 3
        fprintf('[flow] data.%s is %s; collapsing with max over the 3rd dimension\n', ...
            nm, mat2str(size(v)));
        v = max(v, [], 3);
    end
    assert(isequal(size(v), sz), ...
        ['data.%s is %s but the movie is %s. Every mask must match the movie frame -- ' ...
         'a crop applied to one and not the other is the usual cause.'], ...
        nm, mat2str(size(v)), mat2str(sz));
    data.(nm) = v;
end
if strcmp(class(data.bvBW), 'double'), data.bvBW = data.bvBW > 0; end %#ok<STISA>

% ---- per-frame vectors -----------------------------------------------------------------
for f = {'badfr','spike'}
    nm = f{1};
    data = fill(data, nm, false(1,nT));
    v = data.(nm);
    if isempty(v), data.(nm) = false(1,nT); continue, end
    v = v(:)';
    assert(numel(v) == nT, ...
        ['data.%s has %d elements but the movie has %d frames. If the movie was cropped ' ...
         'in time, crop this too.'], nm, numel(v), nT);
    data.(nm) = logical(v);
end

% ---- footprints -------------------------------------------------------------------------
data = fill(data, 'ftprnt', []);
if ~isempty(data.ftprnt)
    assert(isequal(size(data.ftprnt,1), nrow) && isequal(size(data.ftprnt,2), ncol), ...
        'data.ftprnt is %s but the movie frame is %s.', ...
        mat2str(size(data.ftprnt)), mat2str(sz));
    if ~any(data.neuronMap(:))       % derive the projection if it was not supplied
        data.neuronMap = max(double(data.ftprnt), [], 3);
    end
end

% ---- the ANALYSIS MASK ------------------------------------------------------------------
% Every variance figure in the pipeline is computed over these pixels, not the whole
% frame. The frame is mostly background: including it inflates the denominator of every
% "fraction of variance" with pixels that carry no signal, so a regressor that removes a
% lot from the neuron reads as removing a little from the movie. Restricting to the
% footprints makes those fractions mean "of the variance where the neuron actually is".
% PRECEDENCE, highest first:
%   1. data.anaMask, if you set it       -- always wins, nothing is inferred
%   2. any(data.ftprnt > 0, 3)           -- the union of the ROI footprints
%   3. data.neuronMap > 0                -- only sensible if neuronMap IS a footprint
%   4. the whole frame                   -- reported, because the fractions then include
%                                           background and read far too small
data = fill(data, 'anaMask', []);
if isempty(data.anaMask)
    if ~isempty(data.ftprnt)
        data.anaMask = any(double(data.ftprnt) > 0, 3);
    elseif any(data.neuronMap(:))
        % "> 0" is a footprint only if neuronMap IS a footprint. A signed, continuous map
        % -- a spike-triggered image, say -- is positive wherever noise happens to be
        % positive, which is about half the frame. Measured on a synthetic STA: 60% of
        % pixels selected where the neuron covers 21%.
        negFrac = mean(data.neuronMap(:) < 0);
        if negFrac > 0.05
            warning('flow_checkData:signedMask', ...
                ['data.neuronMap has %.0f%% negative pixels, so it looks like a SIGNED map ' ...
                 '(an STA?) rather than a footprint. Thresholding it at 0 selects %.0f%% ' ...
                 'of the frame, which is not a neuron mask. Pass data.ftprnt, or set ' ...
                 'data.anaMask yourself.'], 100*negFrac, 100*mean(data.neuronMap(:) > 0));
        end
        data.anaMask = data.neuronMap > 0;
    else
        data.anaMask = true(sz);
        fprintf(['[flow] no ftprnt and no neuronMap: variance is taken over the whole ' ...
                 'frame, so the fractions include background\n']);
    end
end
assert(isequal(size(data.anaMask), sz), 'data.anaMask is %s but the movie frame is %s.', ...
    mat2str(size(data.anaMask)), mat2str(sz));
data.anaMask = logical(data.anaMask);
assert(any(data.anaMask(:)), 'data.anaMask selects no pixel.');

% ---- traces to project out ---------------------------------------------------------------
data = fill(data, 'voltTr', []);
if ~isempty(data.voltTr)
    V = data.voltTr;
    if size(V,2) ~= nT && size(V,1) == nT
        fprintf('[flow] data.voltTr is [nT x nTrace]; transposed to [%d x %d]\n', ...
            size(V,2), size(V,1));
        V = V.';
    end
    assert(size(V,2) == nT, ...
        ['data.voltTr is %s but the movie has %d frames; it must be [nTrace x nFrames] ' ...
         'or [nFrames x nTrace].'], mat2str(size(data.voltTr)), nT);
    data.voltTr = V;
end
end

function d = fill(d, nm, v)
if ~isfield(d, nm), d.(nm) = v; end
end
