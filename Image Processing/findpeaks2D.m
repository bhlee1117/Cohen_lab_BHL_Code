function [locmask, rois, coids, avgint] = findpeaks2D(img,mindistance,maxpeaks,roiradius)
%findrois   Find 2D peaks and calculate rois around them. 
%   Works well with low noise images of peaks with no background.
%   calculates some statistics of the peaks.
% 
%   2016 Vicente Parot
%   Cohen Lab - Harvard University

    [locmask] = findlocs(img,mindistance,maxpeaks);
    fp8 = find(locmask);
    [~, ix] = sort(img(locmask),'descend');
    sp8 = zeros(size(locmask));
    sp8(fp8(ix)) = 1:numel(fp8);
    dp8 = imdilate(sp8,strel('disk',roiradius));
    rprops = regionprops(dp8,img,'ConvexHull','WeightedCentroid','MeanIntensity');
    rois = {rprops.ConvexHull};
    coids = cell2mat({rprops.WeightedCentroid}');
    avgint = cell2mat({rprops.MeanIntensity}');
end
