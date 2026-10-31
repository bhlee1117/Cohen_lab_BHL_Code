function [w, r] = get_voltPenalty(maps, vmap, mask, lo, hi, wMin)
% [w, r] = get_voltPenalty(maps, vmap, mask, lo, hi, wMin)
%
% Penalise a regressor whose SPATIAL MAP looks like the neuron's voltage footprint.
%
% An oscillation score says a regressor looks like an oscillation; it says nothing about
% WHERE that oscillation lives. A component that is periodic and sits exactly on the
% dendrite is far more likely to be voltage the design has latched onto than blood, and
% regressing it out removes signal. This turns that into a multiplier on the score:
% correlate each regressor's spatial map with the voltage footprint and crush the score
% when they match.
%
%   |r| <= lo          w = 1        (no penalty)
%   lo < |r| < hi      w ramps linearly from 1 down to wMin
%   |r| >= hi          w = wMin     (score effectively destroyed)
%
% The ABSOLUTE correlation is used: a map that is the negative of the footprint is the
% same spatial pattern with the arbitrary sign of the regressor flipped, and is just as
% suspicious.
%
% CHOOSE THE MASK WITH CARE. Correlating over the whole frame is meaningless -- both maps
% are ~0 across the background, which by itself drives r towards 1 for everything. Pass
% the region where the confusion actually arises: the neuron, plus the vessel mask if you
% want vessel pixels to count as evidence against.
%
% A caveat this function cannot resolve: a vessel running ALONG the dendrite overlaps the
% footprint spatially, so a genuine vascular regressor can score a high r. Compare the
% two reference maps ('ftprnt' vs a spike-triggered image) when the penalty fires on a
% regressor you believe is real -- blood is not spike-locked, so the STA image separates
% them where the anatomical footprint cannot.
%
% INPUT
%   maps   [nPix x nReg]  spatial map per regressor (e.g. its covariance with the movie)
%   vmap   [nPix x 1] or [nrow x ncol]   the voltage footprint to compare against
%   mask   logical, [nPix x 1] or [nrow x ncol]   pixels the correlation is taken over
%   lo     |r| below which nothing happens          (default 0.3)
%   hi     |r| at which the full penalty applies    (default 0.6)
%   wMin   the multiplier at full penalty           (default 0.02)
%
% OUTPUT
%   w      [nReg x 1]  score multiplier in [wMin, 1]
%   r      [nReg x 1]  signed correlation with the voltage footprint
%
% 2026.08.28, Byung Hun Lee

if nargin < 4 || isempty(lo),   lo   = 0.3;  end
if nargin < 5 || isempty(hi),   hi   = 0.6;  end
if nargin < 6 || isempty(wMin), wMin = 0.02; end
assert(hi > lo, 'hi (%g) must exceed lo (%g).', hi, lo);

r = get_voltCorr(maps, vmap, mask);       % the correlation itself; this file is the ramp

u = min(max((abs(r) - lo) / (hi - lo), 0), 1);
w = 1 - (1 - wMin) * u;
end
