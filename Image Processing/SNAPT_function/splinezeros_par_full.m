function [delay,F_amp,D_slope,S_width] = splinezeros_par_full(movie, thresh, dir,kernel);
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

% kernel = squeeze(mean(mean(movie, 1), 2));
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
F_amp = zeros(npix, 1);
D_slope = zeros(npix, 1);
int_y = zeros(npix, 1);
parfor j = 1:npix;

    kspline = spapi(3, 1:Lk, kernel);  % spline of the kernel
    % prototype single-pixel functions to be used for evaluation.
    spline1 = kspline;
    dspline1 = fnder(kspline); 

    spline1.coefs = mspline.coefs(j,:);  % copy out the spline coefs for the pixel
    dspline1.coefs = dmspline.coefs(j,:);  % copy out the spline derivative coefs for the pixel
    smax = -fnmin(fncmb(spline1, -1));   % find the maximum of the spline, amplitude
    sbase = mean(fnval(spline1, (Lk-3):Lk));  % average the last 3 frames
    spline2 = fncmb(fncmb(spline1, '-', sbase), 1/(smax - sbase)); % scale the spline
    tmp = fnzeros(fncmb(spline2, '-', thresh));  % find the threshold-crossings
    xsigns = sign(fnval(dspline1, tmp(1,:))); % sort by direction, get value of derivative at crossing point 
    slopeAtFirstCrossing = fnval(dspline1, tmp(1,:));
    %Normalize curve to 1 to get the width (by integration)
    % Find the minimum and substract it to the spline coeff, curve min is
    % at 0
    aux = fncmb(spline1, '-', fnmin(spline1));
    amp_aux = -fnmin(fncmb(aux, -1));   % find the maximum of the spline, again
    normSpline = fncmb(aux,1/amp_aux);
    %integrate
    int_y_aux = integral(@(xx)fnval(normSpline,xx),1,Lk);

    tmp = tmp(1, xsigns*dir > 0);  % only keep the ones in the desired direction
    if ~isempty(tmp)
        [~, indx] = min(abs(tmp - tpeak)); % find the one closest to the peak
        tzero(j) = tmp(indx); %Delay at the zero crossing
        F_amp(j) = smax; % Spline raw amplitude. Delta F
        D_slope(j) = slopeAtFirstCrossing(indx); %Slope at the threshold point
        int_y(j) = int_y_aux;

    else
        tzero(j) =NaN;
        F_amp(j) = NaN;
        D_slope(j) = NaN;
        int_y(j) = NaN;

    end;
%     if round(j/100) == j/100;
%         sprintf(['completed ' num2str(j) ' of ' num2str(npix)])
%     end;
end;
delay = reshape(tzero, [ysize, xsize]);
F_amp = reshape(F_amp, [ysize, xsize]);
D_slope = reshape(D_slope, [ysize, xsize]);
S_width = reshape(int_y, [ysize, xsize]);

% % For testing
% plot( 1:0.1:Lk,fnval(normSpline, 1:0.1:Lk))
% hold on
% plot( 1:0.1:Lk,fnval(dspline1, 1:0.1:Lk))

