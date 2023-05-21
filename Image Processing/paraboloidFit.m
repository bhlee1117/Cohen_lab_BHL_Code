function [out, coeff, center] = paraboloidFit(in)

% [out, coeff, center] = paraboloidFit(in)
%
% This function fits a 2D array to a paraboloid. Useful for image
% background correction.
% Created by JHH 06 March 2009
% Modifed by JHH 20 July 2010 to allow for 'NaN' in image, which is taken
% out
% Modified by JHH 03 August 2010 to calculate and output the center of the
% paraboloid
%
% Details:
% Fits input to equation f(x, y) = A*(x-x0)^2 + B*(y-y0)^2 + C*(x-x0)*(y-y0) + D*(x-x0) +
% E*(y-y0) + F, where (x0, y0) is the center of the input matrix, defined
% as (ceil(sizeX/2), ceil(sizeY/2))
% 
% Input:
% in = 2D array to be fit
%
% Output:
% out = 2D array of size(in), of best fit paraboloid
% coeff = best fit coefficients [A, B, C, D, E, F]
% center = center of paraboloid (column, row)

[sizeY, sizeX] = size(in);

vecX = repmat(floor(-sizeX/2)+1:floor(sizeX/2), sizeY, 1);
vecY = repmat((floor(-sizeY/2)+1:floor(sizeY/2))', 1, sizeX);

basis = [vecX(:).*vecX(:) vecY(:).*vecY(:) vecX(:).*vecY(:) vecX(:) vecY(:) ones(sizeX*sizeY, 1)];

whereNotNan = ~isnan(in);
coeff = basis(whereNotNan, :)\in(whereNotNan); % backslash operator solves system of linear equations for coeff
out = reshape(basis*coeff, sizeY, sizeX); % calculate best fit paraboloid, and reshape to size(in)

% Calculate center of paraboloid by finding where the minimum is located.
% This involves solving the system of equations:
% df/dx = 2A*(x-x0) + C*(y-y0) + D = 0; (1)
% df/dy = 2B*(y-y0) + C*(x-x0) + E = 0; (2)

% Multiply (1) by C and (2) by 2A and subtract:
% C^2*(y-y0) + CD - 4AB*(y-y0) - 2AE = 0
% y = (2AE - CD)/(C^2-4AB) + y0; (3)

% Multiply (1) by 2B and (2) by C and subtract:
% 4AB*(x-x0) + 2BD - C^2*(x-x0) - CE = 0 
% x = (2BD - CE)/(C^2-4AB) + x0; (4)

% Equations 3 and 4 in terms of our variables:
center(1) = (2*coeff(2)*coeff(4) - coeff(3)*coeff(5))/(coeff(3)^2-4*coeff(1)*coeff(2)) + ceil(sizeX/2);
center(2) = (2*coeff(1)*coeff(5) - coeff(3)*coeff(4))/(coeff(3)^2-4*coeff(1)*coeff(2)) + ceil(sizeY/2);





