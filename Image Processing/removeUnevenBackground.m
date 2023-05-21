function [out, background] = removeUnevenBackground(in, border)
% [out, background] = removeUnevenBackground(in, border)
%
% This function corrected for uneven background illumination by fitting a paraboloid and subtracting.
% Created by JHH 20 July 2010
%
% Details:
% Calls paraboloid fit to find best fit paraboloid and subtracts this from
% input image to correct for uneven background illumination.
% Input:
% in = 2D array to be background corrected
% border = number of pixels from the border that the fitting should use. 0
% if whole image is taken.
% 
% Output:
% out = 2D array of size(in), of background corrected image
% background = subtracted paraboloid background;

inOriginal = in;
if border ~= 0
    in(border+1:end-border, border+1:end-border) = NaN;
end

[background, coeff] = paraboloidFit(in);
out = inOriginal - background;