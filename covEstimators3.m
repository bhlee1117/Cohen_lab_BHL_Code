function E = covEstimators3(Xo, Xe, good)
% covEstimators3  Three covariance estimators of the same data, side by side.
%
%   E = covEstimators3(Xo, Xe, good)
%
% Given two simultaneous measurements of the same signal whose noise is
% independent (for voltage imaging: the odd/even checkerboard pixel halves of
% each footprint), returns three estimates of the ROI x ROI covariance:
%
%   (i)   X(t)X(t)        same-time, one channel (the average of the two).
%                         Contains signal + noise. The noise lands ENTIRELY on
%                         the diagonal, inflating every variance and therefore
%                         deflating every correlation.
%   (ii)  Xodd Xeven      same-time, across channels. Independent noise cancels
%                         in expectation, on the diagonal too. Keeps all real
%                         signal at every timescale.
%   (iii) X(t)X(t+1)      one-frame lag, one channel. White noise has zero
%                         lag-1 covariance so this is also unbiased, but it
%                         additionally attenuates real signal that decorrelates
%                         within one frame, by roughly exp(-1/tau).
%
% (ii) and (iii) are both noise-free; (ii) minus (iii) is the real signal living
% at sub-frame timescales. Estimators (ii) and (iii) are NOT positive
% semi-definite: once the signal is exhausted their eigenvalues fluctuate about
% zero and some come out negative. That is correct behaviour for an unbiased
% estimator, not a failure -- E.nPosEig records how many survive.
%
% Inputs
%   Xo, Xe : N x T traces, the two independent-noise channels, same units.
%   good   : 1 x T logical, samples to analyse (blue-on / motion frames out).
%            The lag-1 estimator uses only CONSECUTIVE good pairs.
%
% Output struct E
%   .Cov, .Corr : 1x3 cells, N x N. Corr is Cov normalised by its own diagonal;
%                 rows/columns whose estimated variance is <= 0 become NaN.
%   .Evec, .Eval: 1x3 cells, eigenvectors / eigenvalues sorted descending.
%   .names      : 1x3 labels.
%   .nGood, .nPair, .nROI, .nPosEig, .nBadDiag
%
% See also: cvPCA, get_corrTensor

assert(isequal(size(Xo), size(Xe)), 'covEstimators3:size', 'Xo and Xe must match');
good = logical(good(:)');
assert(numel(good) == size(Xo,2), 'covEstimators3:good', 'good must be 1 x T');

N  = size(Xo,1);
Xc = (Xo + Xe) / 2;                                   % combined channel
g  = good & all(isfinite(Xc),1) & all(isfinite(Xo),1) & all(isfinite(Xe),1);
nG = nnz(g);
assert(nG > N+2, 'covEstimators3:few', 'only %d usable samples for %d ROIs', nG, N);

Cg = Xc(:,g) - mean(Xc(:,g), 2);
Og = Xo(:,g) - mean(Xo(:,g), 2);
Eg = Xe(:,g) - mean(Xe(:,g), 2);

% lag-1 pairs must be CONSECUTIVE in real time, so they are taken on the full
% axis and then restricted to pairs where both frames are analysable
pr = find(g(1:end-1) & g(2:end));
assert(numel(pr) > N+2, 'covEstimators3:fewPairs', 'only %d consecutive good pairs', numel(pr));
Xm = Xc - mean(Xc(:,g), 2);
A1 = Xm(:, pr);  A2 = Xm(:, pr+1);

E.Cov = { Cg*Cg'/(nG-1), ...
          (Og*Eg' + Eg*Og')/2/(nG-1), ...
          (A1*A2' + A2*A1')/2/(numel(pr)-1) };
E.names = {'(i) X(t)X(t) same-time', ...
           '(ii) X_{odd}X_{even} noise-free', ...
           '(iii) X(t)X(t+1) lag-1'};
E.nGood = nG;  E.nPair = numel(pr);  E.nROI = N;

E.Corr = cell(1,3);  E.Evec = cell(1,3);  E.Eval = cell(1,3);
E.nPosEig = zeros(1,3);  E.nBadDiag = zeros(1,3);
for e = 1:3
    C = (E.Cov{e} + E.Cov{e}')/2;                     % enforce exact symmetry
    E.Cov{e} = C;
    d = diag(C);
    E.nBadDiag(e) = nnz(d <= 0);
    sd = sqrt(d);  sd(d <= 0) = NaN;                  % NaN out unusable rows
    E.Corr{e} = C ./ (sd * sd');
    [V, D] = eig(C, 'vector');
    [D, si] = sort(D, 'descend');
    E.Evec{e} = V(:, si);  E.Eval{e} = D;
    E.nPosEig(e) = nnz(D > 0);
end
end
