function [radInfo] = radInterp(img, x0, y0, radii, N)

% [radInfo] = radInterp(img, x0, y0, radii, N)
% Interpolates a 2-d matrix image to find N evenly spaced data values (over 2pi) at fixed radii 
% specified by radii. 
% Created by JHH 03 August 2010, by modifying function radavg
%
% Input:
% img = 2D array to be fit
% x0 = x-coordinate of origin (column)
% y0 = y-coordinate of origin (row)
%
% Output:
% radInfo = structure with three fields:
% .radii: radii at which data points are calculated
% .theta: theta values at which data points are calculated
% .data: 2-d matrix of interpolated data points with 
% rows corresponding to .radii and 
% columns corresponding to .theta

tmpTheta = linspace(0, 2*pi, N+1);
radInfo.theta = tmpTheta(1:end-1);
radInfo.radii = radii';
xi = x0 + radii'*cos(radInfo.theta);
yi = y0 + radii'*sin(radInfo.theta);
radInfo.data = interp2(img, xi, yi, 'cubic');


