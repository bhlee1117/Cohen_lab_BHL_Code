function SI = spatialInfo_mat(SpikeMat, binTr)
% SPATIALINFO_MAT  Vectorized spatial information (bits/spike) for many cells.
%   SI = spatialInfo_mat(SpikeMat, binTr)
%   SpikeMat : nROI x nT  spike counts/rate (NaN allowed)
%   binTr    : 1  x nT  position-bin index per frame; NaN = excluded frame
%              (e.g. non-running).  Must match SpatialInfo.m's convention.
%   SI       : nROI x 1  spatial information (bits per spike).
% This is a batched, allocation-light equivalent of SpatialInfo.m (same result),
% used for the observed value AND every circular-shuffle in S5.
binTr = binTr(:)';
pos_bin = max(binTr);                          % ignores NaN (matches SpatialInfo)
MeanF   = mean(SpikeMat, 2, 'omitnan');        % nROI x 1, over ALL frames
valid   = ~isnan(binTr);                       % frames counted in occupancy denominator
denom   = sum(valid);

b    = binTr(valid);
keep = b>=1 & b<=pos_bin;                       % drop bin 0 / out-of-range (as SpatialInfo does)
cols = b(keep)';
Sv   = SpikeMat(:, valid);
Sv   = Sv(:, keep);                             % nROI x nKeep
Sv(isnan(Sv)) = 0;                              % omitnan means for the per-bin mean

OH     = sparse(1:numel(cols), cols, 1, numel(cols), pos_bin);  % nKeep x pos_bin
binSum = Sv * OH;                                % nROI x pos_bin  (spike sum per bin)
binCnt = full(sum(OH,1));                        % 1 x pos_bin     (frames per bin)

FR = binSum ./ binCnt;                           % nROI x pos_bin  (empty bins -> NaN)
PR = binCnt ./ denom;                            % 1 x pos_bin
SI = sum(PR .* FR .* log2(FR ./ MeanF), 2, 'omitnan') ./ MeanF;   % nROI x 1
end
