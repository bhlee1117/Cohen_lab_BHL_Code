function [Idx, xVal, yVal] = findPoints_mpc(xData, yData, varargin);
% [Idx, xVal, yVal] = findPoints(xData, yData, plotParams);
%
% Plots the xData and yData using usual plotting command syntax, with plot
% parameters specified in strings as in the plot command.
%
% Returns the linear index, xValues, and yValues of points inside a clicked
% region.  The plot should only have one plot on it.  Useful for
% identifying by hand points that have clustered in a scatter plot.
%
% AEC 22 Oct. 2015

%% debugging

% set(0, 'CurrentFigure', figHandle);
% axesObjs = get(figHandle, 'Children');  %axes handles
% dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
% xdata = get(dataObjs(end), 'XData');  %data from low-level grahics objects
% ydata = get(dataObjs(end), 'YData');

plot(xData, yData, varargin{:})
title('Click on the figure to define the ROI')

[xclick, yclick] = getline(gcf);
hold all;
plot([xclick(end); xclick], [yclick(end); yclick]);
hold off;

inside = inpolygon(xData, yData, xclick, yclick);

Idx = find(inside);
xVal = xData(Idx);
yVal = yData(Idx);

