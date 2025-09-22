function [yoffSet, xoffSet] = calculate_shift(avgImg, avgImg2, useManualROI)
% calculate_shift computes the shift to align avgImg with avgImg2
% If useManualROI=true, selects template from avgImg and correlates with avgImg2.
% If useManualROI=false, selects template from avgImg2 and correlates with avgImg.

if nargin < 3
    useManualROI = false;
end

% Normalize
avgImg  = mat2gray(avgImg);
avgImg2 = mat2gray(avgImg2);

if useManualROI
    % User selects ROI from avgImg (template)
    [ROIpoly, ROImask] = get_ROI(avgImg, [], 'Draw ROI from avgImg');
    if isempty(ROImask)
        error('No ROI selected.');
    end
    
    mask = ROImask(:,:,end);
    stats = regionprops(mask, 'BoundingBox');
    bb = round(stats.BoundingBox);  % [x y w h]

    % Crop avgImg and mask to bounding box
    cropped_avgImg = imcrop(avgImg, bb);
    cropped_mask = imcrop(mask, bb);

    % Apply mask after cropping
    template = cropped_avgImg;
    template(~cropped_mask) = 0;

    % Cross-correlate with avgImg2
    cc = normxcorr2(template, avgImg2);

    % Find peak location
    [~, imax] = max(cc(:));
    [ypeak, xpeak] = ind2sub(size(cc), imax);

    % Correct offset calculation
    yoffSet = ypeak - size(template, 1) - bb(2);
    xoffSet = xpeak - size(template, 2) - bb(1);

else
    % Automatic central crop from avgImg2 (template)
    [h2, w2] = size(avgImg2);
    cy = round(h2/2); cx = round(w2/2);
    cropRatio = 0.5;
    halfH = round(h2 * cropRatio / 2);
    halfW = round(w2 * cropRatio / 2);
    template = avgImg2(cy - halfH + 1 : cy + halfH, ...
                       cx - halfW + 1 : cx + halfW);

    bb=[cx - halfW + 1, cy - halfH + 1];
    % Cross-correlate with avgImg
    cc = normxcorr2(template, avgImg);

    % Find peak location
    [~, imax] = max(cc(:));
    [ypeak, xpeak] = ind2sub(size(cc), imax);

    % Correct offset calculation
    yoffSet = -(ypeak - size(template, 1) - bb(2)+2);
    xoffSet = -(xpeak - size(template, 2) - bb(1)+2);


end





end
