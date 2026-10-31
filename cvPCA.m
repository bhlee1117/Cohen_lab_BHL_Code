function [ss, ssRaw, nComp] = cvPCA(X1, X2, nSplit)
% cvPCA  Cross-validated PCA eigenspectrum (Stringer et al. 2019, Fig. 2a).
%
%   [ss, ssRaw, nComp] = cvPCA(X1, X2)
%   [ss, ssRaw, nComp] = cvPCA(X1, X2, nSplit)
%
% Ordinary PCA on a single noisy measurement is BIASED: independent noise adds
% an approximately flat floor to the eigenvalue spectrum, so the tail reports
% noise rather than signal. cvPCA removes that bias by using two measurements
% of the SAME underlying signal whose noise is independent:
%
%       u_n  = n-th left singular vector of X1 on a held-out half of the samples
%       ss_n = < (u_n' * X1) , (u_n' * X2) >   over the OTHER half, matched columns
%
% Because the noise in X_train and X_test is independent it contributes zero in
% expectation to the cross-product, so ss is an unbiased estimate of the signal
% variance in dimension n. Two consequences to expect:
%   * ss_n can be NEGATIVE for large n -- that is the estimator fluctuating
%     around zero once the signal is exhausted, not a bug. Those components are
%     simply where there is no reliable signal left.
%   * the spectrum is only as long as min(size(X,1), rank) -- with a few dozen
%     ROIs you get a few dozen points, far too few to claim a power law (see
%     Stringer et al. Fig. 2g-j, where the usable range grows with unit count).
%
% Inputs
%   X1, X2 : N x T matrices. Same N rows (same ROIs/pixels) and same T samples,
%            measuring the same signal with INDEPENDENT noise. For voltage
%            imaging the natural pair is the odd/even checkerboard pixel halves
%            of each footprint: simultaneous, so any real optical signal
%            (including blood) is shared, while shot noise is independent.
%            Each row is mean-subtracted internally.
%   nSplit : number of random half-splits of the SAMPLES to average over
%            (default 10). Each split also averages both directions
%            (u from X1 applied to X2, and u from X2 applied to X1), which
%            removes the asymmetry of choosing one side as "train".
%
% Outputs
%   ss     : nComp x 1 cross-validated variance per dimension, averaged over
%            splits and normalised to sum to 1 over the POSITIVE components
%            (so it reads as a fraction of reliable variance).
%   ssRaw  : the same before normalisation, in the input units squared.
%   nComp  : numel(ss).
%
% See also: svd, get_eigvector_oddeven

if nargin < 3 || isempty(nSplit), nSplit = 10; end
assert(isequal(size(X1), size(X2)), 'cvPCA:size', 'X1 and X2 must be the same size');

[N, ~] = size(X1);
ok = all(isfinite(X1), 1) & all(isfinite(X2), 1);      % samples usable in both
assert(nnz(ok) > 4, 'cvPCA:noSamples', 'fewer than 5 finite samples in common');
X1 = X1(:, ok);  X2 = X2(:, ok);
Tk = nnz(ok);

X1 = X1 - mean(X1, 2);
X2 = X2 - mean(X2, 2);

nComp = min(N, floor(Tk/2));
acc   = zeros(nComp, 1);
cnt   = 0;

% The cross-product MUST be taken at MATCHING samples: X1(:,t) and X2(:,t)
% observe the same signal at the same instant, and that is the only reason the
% independent noise cancels. Samples are split only to fit u on data that is
% then not evaluated, so u cannot be tuned to the noise it is scored on.
for s = 1:nSplit
    p  = randperm(Tk);
    iFit  = p(1:floor(Tk/2));                          % estimate u here
    iEval = p(floor(Tk/2)+1:end);                      % score here, same columns in both

    for d = 1:2
        if d == 1, F = X1(:, iFit); else, F = X2(:, iFit); end
        F = F - mean(F, 2);
        [U, ~, ~] = svd(F, 'econ');
        k  = min(nComp, size(U, 2));
        E1 = X1(:, iEval) - mean(X1(:, iEval), 2);
        E2 = X2(:, iEval) - mean(X2(:, iEval), 2);
        c1 = U(:, 1:k)' * E1;
        c2 = U(:, 1:k)' * E2;
        v  = sum(c1 .* c2, 2) / (numel(iEval) - 1);    % cross-covariance per dim
        acc(1:k) = acc(1:k) + v;
        cnt = cnt + 1;
    end
end

ssRaw = acc / cnt;
pos   = ssRaw > 0;
if any(pos), ss = ssRaw / sum(ssRaw(pos)); else, ss = ssRaw; end
end
