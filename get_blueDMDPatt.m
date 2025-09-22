function [blueDMDimg, bluePatt, blueCenter]=get_blueDMDPatt(Device_Data,mode)

% [blueDMDimg, bluePatt, blueCenter]=get_blueDMDPatt(Device_Data,mode)
% Description: Generates transformed blue DMD (Digital Micromirror Device) patterns and their centroids 
%              from input device data, applying inverse geometric transformation and cropping to a 
%              specified region of interest (ROI). Supports 'normal' and 'stack' modes for pattern processing.
% Inputs:
%   - Device_Data: Structure containing reference image, transformation data, sensor size, and ROI.
%   - mode: String specifying pattern type ('normal' or 'stack', default: 'normal').
% Outputs:
%   - blueDMDimg: Cropped, transformed DMD pattern image stack.
%   - bluePatt: Cell array of transformed boundary coordinates for each pattern.
%   - blueCenter: Transformed centroid coordinates for each pattern.

if nargin < 2
    mode = 'normal';
end

try
    Rfixed = imref2d(size(Device_Data{1, 6}.refimage.img));
catch
    try
        cam = 2;
        sensorsize = Device_Data{1, 3}.virtualSensorSize; % fusion
    catch
        cam = 1;
        sensorsize = size(Device_Data{1, 3}.testimage); % flash
    end
    Rfixed = imref2d([sensorsize sensorsize]);
end

inverseTform = invert(Device_Data{1, 6}.tform);

switch mode
    case 'normal'
        rawPattern = double(Device_Data{1, 6}.Target);
    case 'stack'
        rawPattern = double(Device_Data{1, 6}.pattern_stack) == 1;
end

% Get centroids from untransformed pattern and transform them
blueCenter = [];
blueDMDimg = [];
bluePatt = cell(size(rawPattern, 3), 1);

for z = 1:size(rawPattern, 3)
    stats = regionprops(rawPattern(:,:,z), 'Centroid');
    boundaries = bwboundaries(rawPattern(:,:,z));
    boundaries= cellfun(@(x) x(:,[2 1]),boundaries,'UniformOutput',false);

    if ~isempty(stats)
        c = stats(1).Centroid;
        blueCenter(z, :) = transformPointsForward(inverseTform, c);
        bluePatt{z} = transformPointsForward(inverseTform, boundaries{1});
    else
        blueCenter(z, :) = [NaN, NaN];
        bluePatt{z}= [NaN, NaN];
    end
end

% Warp the pattern after centroid extraction
revertedImage = imwarp(rawPattern, inverseTform, 'OutputView', Rfixed);

if true
    if cam == 2
        ROI = double(Device_Data{1, 3}.ROI([1 3 2 4]));
    else
        ROI = double(Device_Data{1, 4}.ROI([1 3 2 4]));
    end

    for z = 1:size(revertedImage, 3)
        blueDMDimg(:,:,z) = imcrop(revertedImage(:,:,z), ROI + [0 0 -1 -1]);

        % Convert transformed blueCenter to ROI coordinates
        blueCenter(z, :) = blueCenter(z, :) - [ROI(1)-1, ROI(2)-1];
        bluePatt{z} = bluePatt{z} - [ROI(1)-1, ROI(2)-1];
    end
else
    bluePatt = [];
    blueDMDimg = [];
end
end


