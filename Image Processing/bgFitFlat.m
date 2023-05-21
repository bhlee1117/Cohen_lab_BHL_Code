function [out, coeff] = bgFitFlat(in, bgFxn1, bgFxn2, bgFxn3)

% [out, coeff] = bgFitFlat(in, bgFxn1, bgFxn2, bgFxn3)
% This function fits a 2D array to a linear combination of up to three
% known functions plus a constant flat offset
% Created by JHH 05 April 2011
%
% Details:
% Fits input to equation f(x, y) = A1*bgFxn1(x,y) + A2*bgFxn2(x,y) +
% A3*bgFxn3(x,y) + C
% 
% Input:
% in = 2D array to be fit
% bgFxn1, bgFxn2, bgFxn3 = arrays of size(in)
% (bgFxn2, 
%
% Output:
% out = 2D array of size(in), of best fit f(x,y)
% coeff = best fit coefficients [A1, A2, A3, C]

[sizeY, sizeX] = size(in);

% vecX = repmat(floor(-sizeX/2)+1:floor(sizeX/2), sizeY, 1);
% vecY = repmat((floor(-sizeY/2)+1:floor(sizeY/2))', 1, sizeX);

switch nargin
    case 2
        basis = [bgFxn1(:) ones(sizeX*sizeY, 1)];     
    case 3
        basis = [bgFxn1(:) bgFxn2(:) ones(sizeX*sizeY, 1)];
    case 4
        basis = [bgFxn1(:) bgFxn2(:) bgFxn3(:) ones(sizeX*sizeY, 1)];
end

whereNotNan = ~isnan(in);

coeff = basis(whereNotNan, :)\in(whereNotNan); % backslash operator solves system of linear equations for coeff
% out = reshape(basis*coeff, sizeY, sizeX) + bgDark; % calculate best fit and reshape to size(in)
out = reshape(basis*coeff, sizeY, sizeX); % calculate best fit and reshape to size(in)

