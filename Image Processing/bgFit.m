function [out, coeff] = bgFit(in, bgFxn1, bgFxn2, bgFxn3)

% [out, coeff] = bgFit(in, bgFxn1, bgFxn2, bgFxn3)
% This function fits a 2D array to account for uneven illumination,
% non-parallel flow cells, and dark counts of the camera.
% Created by JHH 06 March 2009
% Modifed by JHH 20 July 2010 to allow for 'NaN' in image, which is taken
% out
% Modified by JHH 04 April 2011 from paraboloidFit
%
% Details:
% Fits input to equation f(x, y) = D*bgFxn1(x,y) + E*bgFxn2(x,y) + A*(x-x0)
% + B*(y-y0) + C
% where (x0, y0) is the center of the input matrix, defined as (ceil(sizeX/2), ceil(sizeY/2))
% 
% Input:
% in = 2D array to be fit
% bgFxn1, bgFxn2, 

% Output:
% out = 2D array of size(in), of best fit paraboloid
% coeff = best fit coefficients [A, B, C, D, E, F]

[sizeY, sizeX] = size(in);

vecX = repmat(floor(-sizeX/2)+1:floor(sizeX/2), sizeY, 1);
vecY = repmat((floor(-sizeY/2)+1:floor(sizeY/2))', 1, sizeX);

switch nargin
    case 2
        basis = [bgFxn1(:) vecX(:) vecY(:) ones(sizeX*sizeY, 1)];     
    case 3
        basis = [bgFxn1(:) bgFxn2(:) vecX(:) vecY(:) ones(sizeX*sizeY, 1)];
    case 4
        basis = [bgFxn1(:) bgFxn2(:) bgFxn3(:) vecX(:) vecY(:) ones(sizeX*sizeY, 1)];
end

whereNotNan = ~isnan(in);

coeff = basis(whereNotNan, :)\in(whereNotNan); % backslash operator solves system of linear equations for coeff
% out = reshape(basis*coeff, sizeY, sizeX) + bgDark; % calculate best fit and reshape to size(in)
out = reshape(basis*coeff, sizeY, sizeX); % calculate best fit and reshape to size(in)

