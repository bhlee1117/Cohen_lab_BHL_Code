function [intensN, pbleach] = rem_pbleach(intens, deterp, prct)

% function [intensN, pbleach] = rem_pbleach_prctile(intens,deterp, prct)
% Removes photobleaching by taking the prct percentile value every deterp and
% smoothing with a window of width deterp.
% This allows removal of negative-going spikes as well as positive-going.
%
% Modified AEC 8 Sept. 2012
% AEC 29 Oct. 2015


intens = intens(:);
L = length(intens);

nPad = deterp - rem(L, deterp);  % Figure out how many elements to pad by
intens = [intens; intens(end)*ones(nPad,1)];
L2 = length(intens);
intens2 = reshape(intens, [deterp, L2/deterp]);

locYVals = prctile(intens2, prct, 1);
locXVals = (1:length(locYVals))*deterp - deterp/2;

pbleach = interp1(locXVals, locYVals, 1:L);

% pbleach = imfilter(pbleach, ones(deterp,1)/deterp, 'replicate');  % This is a way to smooth that handles the boundary conditions better.

% pbleach = smooth(pbleach, deterp);
pbleach = abs(pbleach);
intensN = intens(1:L)./pbleach';


