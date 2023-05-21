function [dfoverfavg, dfoverfbest, baseimg, slopeimg, rsquareimg, outimg] = dfoverf(movie, training, prctilebest, prctilebkg)
% function [dfoverfavg, dfoverfbest, baseimg, slopeimg, rsquareimg, outimg]  = dfoverf(movie, training, prctilebest, prctilebkg)
% extract dFoverF for whole cell and for best pixels, with background
% subtraction
% movie = movie to analyze
% training = signal to use for regression of pixel intensities.  Typically
% training = mean(mean(movie)) .
%
% MAKE SURE TO REMOVE PHOTOBLEACHING FROM MOVIE BEFORE USING THIS FUNCTION.
%
% prctilebest = what top percentile of pixels to use in calculating maximum
% value of dF/F.  Small numbers chooses fewer pixels (more selective).
% prctilebkg = what bottom percentile of pixels to use in calculating the
% background level.  Small numbers chooses fewer pixels (lower background
% value).
% AEC 18 Aug. 2012

avgimg = mean(movie, 3);
[ysize, xsize, nframes] = size(movie);
thresh = prctile(avgimg(:), prctilebkg);
mask = avgimg <= thresh;
bkgintens = mean(mean(mean(repmat(mask, [1 1 nframes]).*movie)))/mean(mask(:));
movie2 = movie - bkgintens;


trainthresh = prctile(training, 5);  % take the bottom 5% of training values as the "zero"
baseindx = find(training < trainthresh);
baseimg = mean(movie2(:,:,baseindx),3);
dFmovie = movie2 - repmat(baseimg, [1 1 nframes]);

dtrain = training - mean(training(baseindx));
dtrain = reshape(dtrain, [1 1 nframes]);
slopeimg = mean(dFmovie.*repmat(dtrain, [ysize, xsize, 1]), 3)/mean(dtrain.*dtrain);  % Linear regression
predmovie = repmat(slopeimg, [1 1 nframes]).*repmat(dtrain, [ysize, xsize, 1]);  % predicted movie
sserr = mean((dFmovie - predmovie).^2, 3);  % Error of the residuals
sstot = mean(dFmovie.^2, 3);  % total variance in data
rsquareimg = 1 - sserr./sstot;

r2threshH = prctile(rsquareimg(:), 100 - prctilebest);
r2maskH = rsquareimg > r2threshH;

dFH = squeeze(mean(mean(dFmovie.*repmat(r2maskH, [1 1 nframes]))))/mean(r2maskH(:));
FH = mean(mean(baseimg.*r2maskH))/mean(r2maskH(:));
dfoverfbest = dFH/FH;
dFL = squeeze(mean(mean(dFmovie)));
FL = mean(baseimg(:));
dfoverfavg = dFL/FL;

outimg = repmat(mat2gray(avgimg), [1 1 3]);
outimg(:,:,2:3) = outimg(:,:,2:3) + repmat(0.2*r2maskH, [1 1 2]);