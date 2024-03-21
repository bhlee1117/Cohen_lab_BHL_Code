function out = splinezeros(movie, thresh, dir);
% function out = splinezeros(movie, thresh, dir);
% uses spline interpolation to find the timing of particular action
% potential features at each pixel in a movie.  The spline-interpolated
% intensity trace at each pixel is mapped so that the max is 1 and the mean
% of the last 10 frames is 0.
% thresh: a number between (0 1)
% Threshold crossings are identified closest to the peak of whole-field
% intensity.
% dir: +1 or -1
% dir sets the direction of the threshold-crossing: +1 = rising edge, -1 =
% falling edge
% AEC 29 Nov. 2013

if isempty(find([-1 1]==dir))
    'dir must be +1 or -1'
    return
end;
if thresh <= 0 || thresh >= 1
    'thresh should be between 0 and 1'
    return
end

kernel = squeeze(mean(mean(movie, 1), 2));
[ysize, xsize, Lk] = size(movie);

mspline = spapi(3, 1:Lk, movie);  % fits a quadratic spline to the data
dmspline = fnder(mspline);  % find the first derivative in time, for assigning directions of zero-crossings
kspline = spapi(3, 1:Lk, kernel);  % spline of the kernel

% prototype single-pixel functions to be used for evaluation.
spline1 = kspline;
dspline1 = fnder(kspline); 

[~, tpeak] = fnmin(fncmb(kspline,-1));  % time of the peak of the kernel
npix = size(mspline.coefs, 1);
tzero = zeros(npix, 1);
for j = 1:npix;
    spline1.coefs = mspline.coefs(j,:);  % copy out the spline coefs for the pixel
    dspline1.coefs = dmspline.coefs(j,:);  % copy out the spline derivative coefs for the pixel
    smax = -fnmin(fncmb(spline1, -1));   % find the maximum of the spline
    sbase = mean(fnval(spline1, (Lk-10):Lk));  % average the last 10 frames
    spline2 = fncmb(fncmb(spline1, '-', sbase), 1/(smax - sbase)); % scale the spline
    tmp = fnzeros(fncmb(spline2, '-', thresh));  % find the threshold-crossings
    xsigns = sign(fnval(dspline1, tmp(1,:))); % sort by direction
    tmp = tmp(1, xsigns*dir > 0);  % only keep the ones in the desired direction
    if ~isempty(tmp)
        [~, indx] = min(abs(tmp - tpeak)); % find the one closest to the peak
        tzero(j) = tmp(indx);
    else
        tzero(j) =NaN;
    end;
    if round(j/100) == j/100;
        sprintf(['completed ' num2str(j) ' of ' num2str(npix)])
    end;
end;
out = reshape(tzero, [ysize, xsize]);

