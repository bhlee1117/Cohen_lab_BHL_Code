function out = spatialfilt(movie, filtsize, sigma);
% function out = spatialfilt(movie, filtsize, sigma);
% returns a spatially filtered movie, where the Gaussian kernel is of dimension
% [filtsize filtsize], with standard deviation sigma.
% if filtsize = 1, then the input movie is returned unmodified.
% The filtering is done with a weighted mean, where each pixel is weighted
% by its standard deviation.  This preferentially weights active pixels in
% the output movie.
% AEC 29 Nov. 2013

stdimg = std(movie, [], 3);
H = fspecial('gaussian', [filtsize, filtsize], sigma);
nframes = size(movie, 3);
out = imfilter(movie.*repmat(stdimg, [1 1 nframes]), H, 'replicate')./repmat(imfilter(stdimg, H, 'replicate'), [1 1 nframes]);