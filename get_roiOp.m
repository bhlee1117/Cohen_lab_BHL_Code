function [Op, nPix] = get_roiOp(ftprnt, ordIdx, mode, thr)
% [Op, nPix] = get_roiOp(ftprnt, ordIdx, mode, thr)
%
% Build the pixels-by-ROI operator that turns a movie into ROI traces:
%
%       T = Op' * tovec(mov)          % nROI x nFrames
%
% so the choice of what "the trace of an ROI" means lives in one place instead of being
% re-derived (and drifted) at every call site.
%
% WHY THE MODE MATTERS. Result.ftprnt is not a 0..1 weight map -- on the 20240807 f=20
% dataset its values run to 8.3e6 with a per-ROI sum of 3.6e8, concentrated in a few
% centre pixels. So the two modes are genuinely different measurements:
%
%   'mean'      plain average over every pixel of the footprint (mask / nPix). Each
%               pixel counts once, so the trace is the ROI's mean dF and dim edge
%               pixels are represented. Noisier per pixel, but not biased toward
%               whatever the brightest pixels happen to contain -- which matters when
%               asking whether a vessel sits on part of an ROI.
%   'weighted'  unit-sum footprint weights, i.e. a weighted average that follows the
%               footprint's own intensity profile. Better SNR, but a handful of centre
%               pixels dominate the result.
%
% Neither is normalised to dF/F; both are linear in the movie, so a linear correction
% applied to the movie can be applied to the trace in closed form.
%
% NOTE ON OVERLAP: footprints are used as they come. Where two ROIs overlap, the shared
% pixels contribute to both traces. Use a label image (max over the 3rd dimension) first
% if disjoint ROIs are needed.
%
% INPUT
%   ftprnt  [H x W x nROIall]  footprints (negatives are clipped to 0)
%   ordIdx  ROI indices to keep, in the order wanted (e.g. Pr.sortDist). [] = all.
%   mode    'mean' (default) or 'weighted'
%   thr     a pixel joins the mask if ftprnt > thr * max(ftprnt) for that ROI.
%           0 (default) = every non-zero pixel. Only used by 'mean'.
%
% OUTPUT
%   Op      [H*W x nROI]  the operator, one column per ROI, columns summing to 1
%   nPix    [nROI x 1]    pixels contributing to each ROI
%
% 2026.08.28, Byung Hun Lee

if nargin < 2 || isempty(ordIdx), ordIdx = 1:size(ftprnt,3); end
if nargin < 3 || isempty(mode),   mode   = 'mean'; end
if nargin < 4 || isempty(thr),    thr    = 0; end

F = double(ftprnt(:,:,ordIdx));
F(F < 0) = 0;
F = reshape(F, [], numel(ordIdx));            % px x nROI

switch lower(mode)
    case 'mean'
        pk = max(F, [], 1);                   % per-ROI peak, for the relative threshold
        M  = F > (thr * pk);                  % logical mask
        nPix = sum(M, 1)';
        assert(all(nPix > 0), ['ROI %s has no pixel above thr = %g of its peak; lower ' ...
            'g7FpThr.'], mat2str(find(nPix == 0)'), thr);
        Op = double(M) ./ nPix';              % each column sums to 1 -> plain mean
    case 'weighted'
        nPix = sum(F > 0, 1)';
        Op   = F ./ max(sum(F, 1), eps);      % each column sums to 1 -> weighted mean
    otherwise
        error('mode must be ''mean'' or ''weighted''.');
end
end
