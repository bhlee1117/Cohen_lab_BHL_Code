function [m, s, n] = binMeanSem(v, bi, nbin)
% BINMEANSEM  Mean, SEM and count of v within integer bins bi (1..nbin).
%   [m, s, n] = binMeanSem(v, bi, nbin)
%   v    : values (any shape)                bi : bin index per value (NaN/out-of-range dropped)
%   nbin : number of bins
%   m,s,n: nbin x 1 mean / standard error of the mean / count per bin (empty bins -> NaN, 0)
    v = v(:); bi = bi(:);
    ok = ~isnan(v) & ~isnan(bi);
    v = v(ok); bi = bi(ok);
    n = accumarray(bi, 1, [nbin 1], @sum, 0);
    m = accumarray(bi, v, [nbin 1], @mean, NaN);
    s = accumarray(bi, v, [nbin 1], @(x) std(x)/sqrt(numel(x)), NaN);
end
