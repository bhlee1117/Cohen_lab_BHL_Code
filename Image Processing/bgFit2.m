function [out, out2, coeff] = bgFit2(in, bgFxn1)

% [out, coeff] = bgFit(in, bgFxn1)
% This function fits a 2D array to account for uneven illumination when it
% is unknown, the theoretical profile of data is known (bgFxn1)
% Created by JHH 06 March 2009
% Modifed by JHH 20 July 2010 to allow for 'NaN' in image, which is taken
% out
% Modified by JHH 04 April 2011 from paraboloidFit
% Modofied bg JHH 11 August 2011 from bgFit

% Details:
% Fits input to equation f(x, y) = (A*(x-x0)^2 + B*(y-y0)^2 + C*(x-x0)*(y-y0) + D*(x-x0) +
% E*(y-y0) + F)*bgFxn1(x,y)
% where (x0, y0) is the center of the input matrix, defined as (ceil(sizeX/2), ceil(sizeY/2))
% 
% Input:
% in = 2D array to be fit
% bgFxn1

% Output:
% out = 2D array of size(in), of best fit
% out2 = 2D array of size(in), of best fit paraboloid that multiplies
% bgFxn1 to get in
% coeff = best fit coefficients [A, B, C, D, E, F]

[sizeY, sizeX] = size(in);

vecX = repmat(floor(-sizeX/2)+1:floor(sizeX/2), sizeY, 1);
vecY = repmat((floor(-sizeY/2)+1:floor(sizeY/2))', 1, sizeX);

basis = [vecX(:).*vecX(:) vecY(:).*vecY(:) vecX(:).*vecY(:) vecX(:) vecY(:) ones(sizeX*sizeY, 1)];
basis2 = basis.*repmat(bgFxn1(:), [1, size(basis, 2)]);
% basis2 = [basis2, ones(sizeX*sizeY, 1)];

whereNotNan = ~isnan(in);
coeff = basis2(whereNotNan, :)\in(whereNotNan); % backslash operator solves system of linear equations for coeff
out = reshape(basis2*coeff, sizeY, sizeX); % 
out2 = reshape(basis*coeff(1:6), sizeY, sizeX); % calculate best fit paraboloid, and reshape to size(in)