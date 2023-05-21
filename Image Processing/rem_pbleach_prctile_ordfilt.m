function [intensN, pbleach] = rem_pbleach_prctile_ordfilt(intensity, npoints, percentile)
% [intensN, pbleach] = rem_pbleach_prctile_ordfilt(intensity, npoints, percentile)
% Removes photobleaching by taking the percentile value every npoints
% with symmetric boundary condition
%
% Modified AEC 8 Sept. 2012
% AEC 29 Oct. 2015
% 2016 Vicente Parot

order = max(1,ceil(npoints*percentile/100));
pbleach = ordfilt2(intensity,order,ones(npoints,1),zeros(npoints,1),'symmetric');
intensN = intensity./pbleach;