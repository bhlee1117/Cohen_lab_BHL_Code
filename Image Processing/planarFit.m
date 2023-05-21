function [out, coeff] = planarFit(in)

% [out, coeff] = planarFit(in)
%
% This function fits a 2D array to a plane. Useful for image
% background correction.
% Created by JHH 06 March 2009
% Modified by JHH 20 July 2010 to allow for 'NaN' in image, which is taken
% out
% Modified by JHH 22 March 2011 from paraboloidFit for planarFit
%
% Details:
% Fits input to equation f(x, y) = D*(x-x0) + E*(y-y0) + F, where (x0, y0) is the center of the input matrix, defined
% as (ceil(sizeX/2), ceil(sizeY/2))
% 
% Input:
% in = 2D array to be fit
%
% Output:
% out = 2D array of size(in), of best fit paraboloid
% coeff = best fit coefficients [A, B, C, D, E, F]

[sizeY, sizeX] = size(in);

vecX = repmat(floor(-sizeX/2)+1:floor(sizeX/2), sizeY, 1);
vecY = repmat((floor(-sizeY/2)+1:floor(sizeY/2))', 1, sizeX);

basis = [vecX(:) vecY(:) ones(sizeX*sizeY, 1)];

whereNotNan = ~isnan(in);
coeff = basis(whereNotNan, :)\in(whereNotNan); % backslash operator solves system of linear equations for coeff
out = reshape(basis*coeff, sizeY, sizeX); % calculate best fit paraboloid, and reshape to size(in)
