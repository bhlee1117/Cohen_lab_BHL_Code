function [Idx, xVal, yVal] = nested_findPoints(xData, yData, varargin);
% [Idx, xVal, yVal] = nested_findPoints(xData, yData, plotParams);
%
% Plots the xData and yData using usual plotting command syntax, with plot
% parameters specified in strings as in the plot command.
%
% Returns the linear index, xValues, and yValues of points inside a clicked
% region.  The plot should only have one plot on it.  Useful for
% identifying by hand points that have clustered in a scatter plot.
% Modified by Tian from findPoints
% Tian He 02/14/2019


 %% 
plot(xData, yData, varargin{:});
title('Click on the figure to zoom in')
ax = gca;
xlimits=ax.XLim;
ylimits=ax.YLim;
rois_points={};
npts=1;
counter=1;
Idx=[];
xVal=[];
yVal=[];
 while(npts > 0)
    rectpos = getrect;
    xv = rectpos(ones(5,1),1) + [0 rectpos([3 3]) 0 0]';
    yv = rectpos(ones(5,1),2) + [0 0 rectpos([4 4]) 0]';
    if size(unique([xv, yv],'rows'),1) == 1, break, end
    in = inpolygon(xData, yData,xv,yv);

    plot(xData(find(in==1)), yData(find(in==1)), varargin{:});
    drawnow


    [xclick, yclick] = getline(gcf);
    rois_points{counter}=[xclick, yclick];
    counter=counter+1;


    plot(xData, yData, varargin{:});
    hold on
    plot([xclick(end); xclick], [yclick(end); yclick]);
    hold off;
    drawnow

    inside = inpolygon(xData, yData, xclick, yclick);

    Idx = [Idx,find(inside)'];
    xVal =[xVal,xData(Idx)'];
    yVal =[yVal, yData(Idx)'];
 end
Idx=unique(Idx);
xVal=unique(xVal);
yVal=unique(yVal);

clf;
plot(xData, yData, varargin{:});
hold on
for i=1:numel(rois_points)
        temp=rois_points{i};
        clear xclick y click;
        xclick=temp(:,1);
        yclick=temp(:,2);
        plot ([xclick(end); xclick], [yclick(end); yclick],'-r','LineWidth',1);
        hold on
 end
hold off
axis([xlimits,ylimits]);





