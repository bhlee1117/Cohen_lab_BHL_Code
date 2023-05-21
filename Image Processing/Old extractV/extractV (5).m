function [Vout, corrimg, weightimg, offsetimg] = extractV(imgs, Vin)
% [Vout, corrimg, weightimg, offsetimg] = extractV(imgs, Vin)
% 
% Identifies and weights voltage-responsive pixels in a movie of a cell
% expressing a voltage indicating fluorescent protein.  Constructs an
% estimate of the membrane potential from the pixel weight matrix.
% 
% Inputs
% imgs: [M N K] array of images (size M by N), with K frames in the movie.
% Vin: K by 1 vector containing the measured membrane potential at each
% frame.  In the absence of electrical recording, Vin can be replaced by
% the whole-cell fluorescence.  Running extractV iteratively then
% identifies pixels whose intensity co-varies with the mean.  This
% procedure is useful for spike identification, but the offset and scale of
% Vout are then arbitrary.
% 
% extractV uses a linear model of the fluorescence vs. voltage at each
% pixel, i:
% S(i,t) = V(t)*m(i) + b(i) + e(i,t)
% where S(i,t) is the movie, m(i) is a position-dependent weight factor,
% b(i) is a position-dependent offset and e(i,t) is the noise, assumed to
% be spatially and temporally uncorrelated Gaussian white noise, with a
% pixel-specific variance.
%
% Outputs
% Vout: K by 1 vector of estimated membrane voltages
% corrimg: Image showing correlation between Vin and signal at each pixel.
% weightimg: Weighting coefficients assigned to each pixel, based on
% correlation with Vin and residual noise.
% offsetimg: Offset to be added to each pixel.
%
% Caveats: 
% 1) Assumes voltage indicator responds instantaneously to changes in
% V_m.  This is not a temporal deconvolution algorithm.
% 2) Assumes no time delays in voltage propagation throughout or between
% cells.
% 3) Does not handle motion artifacts or photobleaching.  If either of
% these is a possibility, the data should first be corrected by
% translating, morphing, and normalizing to a photobleaching fit.
% 4) Beware of over fitting!  If the movie has a lot of noise, and the
% number of pixels (M x N) is large compared to the number of frames (K),
% then there will exist a linear combination of pixel intensities that
% yields a Vout that faithfully reproduces Vin.  However, the weightimg
% will look like noise.  So it is important to check that the weightimg is
% consistent with expectations.
%
% Adam E. Cohen 22 Jan 2011

avgimg = mean(imgs, 3);
avgV = mean(Vin);
dV = Vin - avgV;
L = length(Vin);

% subtract off background
imgs = imgs - repmat(avgimg, [1 1 L]);

% correlate the changes in intensity with the applied voltage
[ysize, xsize] = size(avgimg);

dV2 = reshape(dV, 1, 1, L);
dVmat = repmat(dV2, [ysize xsize 1]);
corrimg = mean(dVmat.*imgs,3);
corrimg = corrimg/mean(dV.^2);

% calculate a dV estimate at each pixel, based on the linear regression.
imgs2 = zeros(size(imgs));
corrmat = repmat(corrimg, [1 1 L]);
imgs2 = imgs./corrmat;
clear corrmat

% Look at the residuals to get a noise at each pixel
sigmaimg = mean((imgs2 - dVmat).^2,3);

weightimg = (1./sigmaimg);
weightimg = weightimg/mean(weightimg(:));

dVout = squeeze(mean(mean(imgs2.*repmat(weightimg, [1 1 L]))));

Vout = dVout + avgV;
offsetimg = avgimg - avgV*corrimg;

    