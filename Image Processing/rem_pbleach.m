function [intensN, pbleach] = rem_pbleach(intens, deterp)

% function [intensN, pbleach] = rem_pbleach(intens,deterp)
% Removes photobleaching by taking the minimum value every deterp and
% smoothing with a window of width deterp
% Modified AEC 8 Sept. 2012


intens = intens(:);

pbleach = imerode(intens, ones(deterp,1));
pbleach = imfilter(pbleach, ones(deterp,1)/deterp, 'replicate');  % This is a way to smooth that handles the boundary conditions better.

% pbleach = smooth(pbleach, deterp);
pbleach = abs(pbleach);
intensN = intens./pbleach;


% nframes = length(intens);
% ndeterp = nframes/deterp;
% tmpD = zeros(ndeterp,1); timeD = zeros(ndeterp,1);
% for j = 1:ndeterp
%     indx = (j-1)*deterp+1:j*deterp;
%     tmpD(j) = min(intens(indx));
%     timeD(j) = find(intens == tmpD(j),1,'last');
% end
% pbleach = interp1(timeD,smooth(tmpD,5),1:nframes,'linear','extrap');
% intensN = intens./pbleach';

