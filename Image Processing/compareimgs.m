function [out, corrimg] = compareimgs(A, B, baseline);
% function [out, corrimg] = compareimgs(A, B, baseline);
% Compares two images, A and B and returns a similarity score.  Allows for
% translation between the images.  Useful for identifying images of the
% same cell from different times.
% AEC 17 Dec. 2011

Aavg = mean(A(:));
Bavg = mean(B(:));

A = A - Aavg;
B = B - Bavg;

% Aflat = ones(size(A));
% Bflat = ones(size(B));
% 
% baseline = xcorr2(Aflat, Bflat);

ABcorr = conv2(A, fliplr(flipud(B)), 'same');
ABcorr = ABcorr./baseline;
corrimg = ABcorr/(std(A(:))*std(B(:)));
out = max(corrimg(:));

