function [intensN, pbleach] = rem_pbleach_iterate(intens, deterp, nIteration)

% function [intensN, pbleach] = rem_pbleach_iterate(intens,deterp,nIteration)
% Removes photobleaching by taking the minimum value every deterp and
% smoothing with a window of width deterp.
%
% Because it doesn't deal well with steep traces, expecially when using a
% large deterp window, repeat nIteration times to iteratively flatten the baseline.
% Modified AEC 8 Sept. 2012
% KAW Mar. 2016

intens = intens(:);     % convert to a column vector

pbleach = zeros(size(intens));
intensN = intens;
for i = 1:nIteration
    bkgd = imerode(intensN, ones(deterp,1));
    bkgd = imfilter(bkgd, ones(deterp,1)/deterp, 'replicate');  % This is a way to smooth that handles the boundary conditions better.
    bkgd = bkgd - min(bkgd);
    intensN = intensN - bkgd;
    pbleach = pbleach + bkgd;
end
pbleach = pbleach + min(intensN);   % add the offset