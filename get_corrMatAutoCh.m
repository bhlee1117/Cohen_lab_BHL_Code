function C = get_corrMatAutoCh(sub_comb, sub_ch1, sub_ch2, frameRange)
% GET_CORRMATAUTOCH  Pairwise ROI correlation matrix with a noise-free diagonal.
%
%   C = get_corrMatAutoCh(sub_comb, sub_ch1, sub_ch2, frameRange)
%
% Correlations between different ROIs (off-diagonal) are computed from the
% combined trace SUB_COMB. The auto-correlation of each ROI (diagonal) is
% computed by correlating the two independent optical channels SUB_CH1 and
% SUB_CH2: their shot noise is independent, so the zero-lag value is not
% inflated by noise the way a single-trace autocorrelation would be.
%
% Inputs:
%   sub_comb          - N x T combined subthreshold trace (NaN allowed)
%   sub_ch1, sub_ch2  - N x T per-channel subthreshold traces (NaN allowed)
%   frameRange        - columns (frames) to include
%
% Output:
%   C                 - N x N correlation matrix

C = get_corrMat(sub_comb, sub_comb, frameRange);   % off-diagonal (different ROIs)
A = get_corrMat(sub_ch1, sub_ch2, frameRange);     % two-channel auto (same ROI)
d = logical(eye(size(C)));
C(d) = A(d);
end
