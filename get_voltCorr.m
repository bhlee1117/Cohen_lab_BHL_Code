function r = get_voltCorr(maps, vmap, mask)
% r = get_voltCorr(maps, vmap, mask)
%
% Correlate each regressor's SPATIAL MAP with a voltage footprint, over a chosen mask.
%
% This is the measurement on its own, with no policy attached. What a caller does with a
% high |r| -- ramp a penalty (get_voltPenalty) or veto the row outright (flow_regressors)
% -- is the caller's business.
%
% The correlation is signed, but callers should judge |r|: a map that is the negative of
% the footprint is the same spatial pattern with the arbitrary sign of the regressor
% flipped, and is just as suspicious.
%
% CHOOSE THE MASK WITH CARE. Correlating over the whole frame is meaningless -- both maps
% are ~0 across the background, which by itself drives r towards 1 for everything. Pass
% the region where the confusion actually arises: the neuron, plus the vessel mask if you
% want vessel pixels to count as evidence against.
%
% INPUT
%   maps   [nPix x nReg]  spatial map per regressor (e.g. its covariance with the movie)
%   vmap   [nPix x 1] or [nrow x ncol]   the voltage footprint to compare against
%   mask   logical, [nPix x 1] or [nrow x ncol]   pixels the correlation is taken over
%
% OUTPUT
%   r      [nReg x 1]  signed correlation with the voltage footprint
%
% 2026.08.30, Byung Hun Lee

idx = find(mask(:));
assert(~isempty(idx), 'The mask selects no pixel.');

M = double(maps(idx,:));
v = double(vmap(:));  v = v(idx);

M = M - mean(M, 1);
v = v - mean(v);
den = sqrt(sum(M.^2, 1))' * sqrt(sum(v.^2));
r   = (M' * v) ./ max(den, eps);
r(~isfinite(r)) = 0;                      % a constant map correlates with nothing
end
