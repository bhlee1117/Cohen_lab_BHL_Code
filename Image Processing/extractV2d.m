function [Vout, corrimg, weightimg, offsetimg] = extractV2d(imgs, Vin)
% [Vout, corrimg, weightimg, offsetimg] = extractV2d(imgs, Vin)
% 
% Identifies and weights voltage-responsive pixels in a movie of a cell
% expressing a voltage indicating fluorescent protein.  Constructs an
% estimate of the membrane potential from a weighted sum of the intensities
% of each pixel in each frame of the movie.
% 
% Inputs
% imgs: [M K] array of images (size M by 1), with K frames in the movie.
% Vin: K by 1 vector containing the measured membrane potential at each
% frame.  In the absence of electrical recording, Vin can be replaced by
% the whole-cell fluorescence.  Running extractV iteratively then
% identifies pixels whose intensity co-varies with the mean.  This
% procedure produces a weight matrix nearly identical to the one produced 
% with a real patch recording, which is useful for spike identification.
% However, if extractV is trained with whole-field intensity, then the 
% offset and scale of Vout are arbitrary.
% 
% extractV uses a linear model of the fluorescence vs. voltage at each
% pixel, i:
% S(i,t) = V(t)*m(i) + b(i) + e(i,t)
% where S(i,t) is the movie, m(i) is a position-dependent sensitivity,
% b(i) is a position-dependent offset and e(i,t) is the noise, assumed to
% be spatially and temporally uncorrelated Gaussian white noise, with a
% pixel-specific variance.
%
% Outputs
% Vout: K by 1 vector of estimated membrane voltages.  Vout is the
% least-square estimate, assuming the model above.
% corrimg: Image showing correlation between Vin and signal at each pixel.
% This image alone is insufficient to determine the weight matrix because
% it does not take into account the level of noise at each pixel.
% weightimg: Weighting coefficients assigned to each pixel, based on
% correlation with Vin and residual noise.
% offsetimg: Offset to be added to each pixel to produce a Vout with the
% correct offset.
%
% Caveats: 
% 1) Assumes voltage indicator responds instantaneously to changes in
% V_m.  This is not a temporal deconvolution algorithm.
% 2) Assumes no time delays in voltage propagation throughout or between
% cells.
% 3) Does not handle motion artifacts or photobleaching.  If either of
% these is a possibility, the data should first be corrected by
% translating, morphing, and normalizing to a photobleaching fit.
% 4) Assumes only a single value for the membrane potential.  If different
% cells are doing different things, the image should be segmented first, or
% a more sophisticated algorithm such as Independent Component Analysis
% should be used.
% 4) Beware of over fitting!  If the movie has a lot of noise, and the
% number of pixels (M x 1) is large compared to the number of frames (K),
% then there will exist a linear combination of pixel intensities that
% yields a Vout that faithfully reproduces Vin.  However, the weightimg
% will look like noise.  It is important to check that the weightimg is
% consistent with expectations.  Alternatively one may train extractV on
% one piece of a data set, and then apply it to a different piece, using
% ApplyWeights.m
%
% Adam E. Cohen 22 Jan 2011
% JHH 28 March 2013, edited to accomodate 1d image, 2d data array of MxK.

avgimg = mean(imgs, 2);
avgV = mean(Vin);
dV = Vin - avgV;
L = length(Vin);

% subtract off background
imgs = imgs - repmat(avgimg, [1 L]);

% correlate the changes in intensity with the applied voltage
[ysize, xsize] = size(avgimg);

dV2 = reshape(dV, 1, L);
dVmat = repmat(dV2, [ysize xsize]);
corrimg = mean(dVmat.*imgs,2);
corrimg = corrimg/mean(dV.^2);

% calculate a dV estimate at each pixel, based on the linear regression.
imgs2 = zeros(size(imgs));
corrmat = repmat(corrimg, [1 L]);
imgs2 = imgs./corrmat;  % This will give NaN where corrmat == 0.
clear corrmat

% Look at the residuals to get a noise at each pixel
sigmaimg = mean((imgs2 - dVmat).^2,2);

weightimg = (1./sigmaimg);
weightimg(isnan(weightimg(:))) = 0;
weightimg = weightimg/mean(weightimg(:));
imgs2(isnan(imgs2(:))) = 0;  % Set places where imgs2 == NaN to zero
dVout = mean(imgs2.*repmat(weightimg, [1 L]));

Vout = dVout + avgV;
offsetimg = avgimg - avgV*corrimg;

    